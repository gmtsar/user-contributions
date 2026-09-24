# xcorr3

[English](README.md)

如果你已经在用 GMTSAR 的 `xcorr`，可以用 xcorr3 选择 CPU 多进程计算（`mt`）
或 NVIDIA GPU 计算（`cc`）。输入仍是两份 PRM 文件，采样、搜索和插值参数也
沿用原来的调用方式。需要对比结果时，可以切回原版 `xcorr`；不指定时默认使用原版。

## 安装

下面的命令在 Linux 中执行，WSL2 也可以。先准备好 Python 3 和可用的 GMTSAR。
CPU 版本还需要 GCC/Make，以及与 GMTSAR 构建配套的 GMT、LAPACK、BLAS 和
libtiff 开发库。

在本目录执行，把 `GMTSAR_HOME` 改成你已构建的 GMTSAR 源码目录：

```sh
make install-mt GMTSAR_HOME=/usr/local/GMTSAR PREFIX="$HOME/.local"
export PATH="$HOME/.local/bin:$PATH"
```

这样就装好了 `xcorr3` 和 `xcorr_mt`。把 PATH 那一行加入 shell 启动文件，
以后打开终端就能直接使用。原有的 `xcorr` 会保留。

如果要用 GPU，还需 NVIDIA 驱动、CUDA Toolkit、兼容的 C/C++ 编译器、Make、
pkg-config 和 GLib 开发库，然后执行：

```sh
make cc                              # CUDA 位于 /usr/local/cuda
make install-cc PREFIX="$HOME/.local"
```

如果用的是 `/usr/bin` 下的系统 CUDA，将 `make cc` 换成 `make cc-system`。
只安装入口则执行 `make install PREFIX="$HOME/.local"`。

## 试着运行

进入数据目录，把下面的文件名换成你的两份 PRM。PRM 中记录的 SLC 路径需要
能从当前目录访问。

```sh
xcorr3 master.PRM secondary.PRM                          # GMTSAR 原版
xcorr3 master.PRM secondary.PRM --backend mt -nproc 6     # CPU，6 个进程
xcorr3 master.PRM secondary.PRM --backend cc             # NVIDIA GPU
```

需要调整采样密度，可以加上 `-nx 20 -ny 50`，分别在距离向和方位向取 20、50 个
采样位置；搜索半径用 `-xsearch 128 -ysearch 128` 设置。更多参数见 `xcorr3 --help`，
少传输入文件时也会显示用法。

通常会得到 `freq_xcorr.dat`：原版和 MT 为五列，CC 多一列 `peak_snr`。
CC 用于 SLC 的频域互相关，结果可能与 CPU 版本不同。我们用定日地震数据做了
对比，[验证说明中附有结果和图件](VALIDATION.md)。

如果想在已有 GMTSAR 配置中选择计算方式，或接入处理脚本，见[详细用法](USAGE.zh-CN.md)。
GMTSAR 的 TOPS 几何配准和 ESD 不调用 `xcorr`，因此不会由这个入口加速。

入口和 MT 使用 [GPL-3.0-or-later](LICENSE)，MT 基于 GMTSAR `xcorr`；
CC 保留其 [MIT 许可及 CUI Hao/Jazz-0626 版权声明](backends/cc/LICENSE)。
