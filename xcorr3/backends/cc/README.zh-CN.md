# xcorr_cc

[English](README.md)

基于 CUDA/cuFFT 的 SAR 振幅互相关程序，用于影像配准与像素偏移追踪。
读取两组 GMTSAR PRM/SLC，输出采样点的距离向、方位向偏移；可选地理编码。

## 依赖与安装

- Linux、NVIDIA GPU 及兼容驱动、包含 `nvcc`/cuFFT 的 CUDA Toolkit。
- Make、与 CUDA 兼容的 C/C++ 编译器、`pkg-config`、GLib 开发库。
- 仅 `-geocode` 额外需要 GMT、GMTSAR（`proj_ra2ll.csh` 在 `PATH` 中）、
  C shell，以及 GMTSAR 生成的匹配 `trans.dat`；仅互相关不需要这些依赖。

在源码目录执行（Ubuntu/Debian 示例，CUDA 需另行安装）：

```sh
sudo apt install build-essential pkg-config libglib2.0-dev
make                         # CUDA 位于 /usr/local/cuda
# 若使用 /usr/bin 下的系统 CUDA，改用：bash build_local.sh
mkdir -p "$HOME/.local/bin"
install -m 755 xcorr_cc "$HOME/.local/bin/xcorr_cc"
export PATH="$HOME/.local/bin:$PATH"   # 写入 shell 启动文件可永久生效
```

可手动指定架构，例如 RTX 2070 SUPER 使用 `make CUDA_ARCH=sm_75`。
其他 CUDA 路径可设置 `CUDA_HOME` 或 Makefile 的 `CUDA_NVCC`、`CUDA_INC`、
`CUDA_LIB`。本地回归环境：CUDA 12.0、RTX 2070 SUPER、WSL2 Ubuntu 24.04。

## 使用示例

在两份 PRM 所记录的 SLC 路径均可访问的目录中执行：

```sh
# 仅互相关，保留 PRM 中的粗偏移
xcorr_cc master.PRM secondary.PRM -nx 20 -ny 50 -xsearch 128 -ysearch 128

# 已配准影像对的可选地理编码；当前目录需有匹配的 trans.dat
xcorr_cc master.PRM aligned.PRM -nx 40 -ny 100 -xsearch 128 -ysearch 128 \
  -noshift -geocode -snr 18 -psnr 5

xcorr_cc --help
```

示例参数需按数据调整。`-noshift` 会忽略 PRM 粗偏移，应按输入状态明确选择。

| 参数 | 默认值 | 用途 |
|---|---|---|
| `-nx N`、`-ny N` | 16、32 | 距离向、方位向采样点数 |
| `-xsearch N`、`-ysearch N` | 64、64 | 搜索半径，须为 2 的幂；启用亚像素插值时至少为 8 |
| `-range_interp N`、`-interp N` | 2、16 | 距离向、亚像素插值倍数；距离向倍数须为 2 的幂 |
| `-norange`、`-nointerp` | 关闭 | 禁用对应插值 |
| `-noshift` | 关闭 | 忽略 PRM 的 `rshift`/`ashift` |
| `-geocode` | 关闭 | 过滤、网格化、转米及地理编码 |
| `-snr N`、`-psnr N` | 10、5 | 地理编码的相关系数、峰显著性阈值 |
| `-no_blockmedian` | 关闭 | 直接网格化，不做中值处理 |

`freq_xcorr.dat` 为六列：
`x_pixel x_offset y_pixel y_offset correlation peak_snr`。
偏移单位为像素，相关系数范围为 0–100。`-geocode` 额外输出
`azi_offset.grd`、`rng_offset.grd` 及其 `*_offset_ll.grd` 地理坐标网格，单位为米；
集成后处理的距离向转换使用像素偏移的**负值**。

本程序不覆盖 GMTSAR `xcorr` 的全部模式，不支持 CPU 时域与网格输入模式，
CUDA 与 CPU 结果可能存在差异。请检查退出码：输入、输出或请求的地理编码失败
均返回非零。结果先暂存，检测到发布失败会回滚；同一输出目录只运行一个任务，
多个文件的发布不保证强制终止时整体原子性。测试方法见 [测试说明](test/README.md)。

## 许可与致谢

[MIT 许可](LICENSE)。基于 CUI Hao 的 *Parallel xcorr programs for GMTSAR*
及 GMTSAR 互相关算法；CUDA/cuFFT 移植与 POT 脚本由 Jazz-0626 提供。
版权信息见 `LICENSE`。
