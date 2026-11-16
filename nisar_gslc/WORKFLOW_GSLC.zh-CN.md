# 处理流程与参数

[English](WORKFLOW_GSLC.md) | **简体中文** · [首页](README.zh-CN.md)

以下命令均在 `nisar_gslc/` 包目录内运行。

## 1. 准备项目

在 `raw/` 放入两景包含 A/B 的 GSLC 后执行：

```bash
./run1_prepare_GSLC_LA_LB.csh
```

准备步骤读取成像时间、频率、投影和像元间距，检查数据对，建立 LA/LB、配置和
`gslc_project.json`。原始 HDF5 只读。默认处理整幅共同地图网格，具体有效覆盖由产品掩膜确定。

小范围诊断可以在准备时加 `--bbox XMIN XMAX YMIN YMAX`，坐标为产品 CRS 的米。
换数据对或范围应新建项目。配对名称使用 UTC 日期和标准年内日序（1 月 1 日为 001），
可能与旧 GMTSAR 零起算日序的目录名称不同。

## 2. 修改配置

准备后实际使用：

```text
LA/config_GSLC_A.txt
LB/config_GSLC_B.txt
config_IONO.txt
```

格式为 `参数名 = 数值`，`#` 后为注释。`configs/` 模板只在新项目准备时复制；
修改模板不会修改已有项目的实际配置。

| 参数 | A 默认 | B 默认 | 含义 |
| --- | --- | --- | --- |
| `filter_wavelength` | 120 | 480 | 干涉高斯滤波波长，m |
| `output_spacing_x` | 40 | 120 | 东西向网格间距，m |
| `output_spacing_y` | 40 | 120 | 南北向网格间距，m |
| `goldstein_mode` | coherence | coherence | 按窗口相干性自适应指数 |
| `filter_patch` | 32 | 32 | Goldstein 窗口边长，输出像元 |
| `correlation_threshold` | 0.1 | 0.1 | 解缠质量阈值 |
| `unwrap_mode` | smooth | smooth | SNAPHU 代价模式 |
| `effective_looks` | 32 | 8 | 暂定的 SNAPHU 统计参数 |
| `unwrap_tiles_y/x` | 4 / 4 | 2 / 2 | 整幅分块，南北/东西 |
| `unwrap_overlap_y/x` | 512 / 512 | 512 / 512 | 分块重叠，输出像元 |
| `unwrap_nproc` | 2 | 2 | 同一频段内并行求解块数 |
| `min_valid_fraction` | 0.8 | 0.8 | 原始块/高斯权重的最低有效比例 |

输出间距必须分别是 HDF5 原生间距的整数倍，不会将 B 上采样成不存在的更细信息。
高斯波长独立于网格间距；标准差沿用 GMTSAR 的 `sigma = 波长 / 5.3`，支持范围约为 ±2 sigma。

固定 Goldstein 指数时设 `goldstein_mode = fixed` 和 `filter_alpha = 0.5`。
自适应模式保持 `filter_alpha` 为注释。`filter_patch` 应是至少为 8 的 2 次幂。
小诊断窗口可能需要在两频段配置中将 `unwrap_tiles_y/x` 均改成 1。
A/B 顺序处理，`unwrap_nproc` 控制每个频段内部的分块并行。

独立的电离层配置为：

```text
filter_wavelength = 20000
prefilter_a = 1000
coherence_threshold = 0.1
coarse_samples_per_wavelength = 8
surface_tension = 0.1
interpolate_phase = 0
nearest_grid_radius = 300
```

这里 `filter_wavelength` 和 `prefilter_a` 的单位均为 m，与 A/B 干涉高斯滤波分开。
`interpolate_phase = 1` 才启用 run4/run6 的求解器输入最近邻补相位；半径单位是输出像元。
补值仅辅助求解，正式结果恢复原始有效性/相干性掩膜。默认 0，不补相位。

## 3. 干涉与解缠

Run2 根据实际地图像元中心匹配两景，核对 CRS、间距和网格兼容性，读取原生有效性掩膜，
只使用两景共同有效的复数像元。相位约定为：

```text
I = reference * conjugate(secondary)
P1 = abs(reference)**2
P2 = abs(secondary)**2
phase = angle(weighted_mean(I))
coherence = abs(weighted_mean(I)) / sqrt(weighted_mean(P1) * weighted_mean(P2))
```

高斯权重按有效支持比例归一化，再按目标间距采样。相干性随后使用类似 GMTSAR 的 3×3 权重，
原有掩膜空洞继续保留。`phase.grd` 为高斯加权后的包裹相位，`corr.grd` 为相干性。
调用已有 `phasefilt` 完成 Goldstein，生成 `phasefilt.grd`；自适应指数为
`1 - patch_weighted_coherence`。此前仅用于对照的 `phase_raw.grd` 不在分步流程中生成。

Run3 预览质量掩膜。Run4 读取已有滤波相位和相干性，运行 SNAPHU，恢复通过筛选的有效像元，
输出 `unwrap.grd`、`conncomp.grd`、`los.grd`。默认不生成陆海掩膜。
如果频段结果目录里有可选 `landmask.grd`，必须与结果完全同网格，大于 0 的位置为保留。

Run4 复用完整且未过时的 `display_amp.pdf`、`phase.pdf`、`corr.pdf`、`phasefilt.pdf` 及 CPT，
新增 `unwrap.pdf`、`los.pdf`，更新完整图集；缺图或过时的图自动补画。
`plots.json` 的 `reused_pdfs`、`generated_pdfs` 记录两类图。需要强制重画时：

```bash
python plot_gslc_gmt.py LA/intf/<pair> --overwrite
```

将 `<pair>` 替换为实际配对目录。

## 4. 电离层估计与 A 再次解缠

Run5 按 B 的实际地图坐标线性插值 A 解缠相位，最近邻采样 A 相干性，并按配置对 A 在 B 网格上预滤波。
频率从频段报告记录的 HDF5 元数据读取。`fB > fA` 时，空间滤波前的 A 频率电离层估计为：

```text
r = fA / fB
iono_A = (phase_A_on_B - r * phase_B) / (1 - r*r)
```

按 A/B 平均相干性和频率相关离群区间选择可靠样本，沿用原方法进行最近邻填补，以及
median → boxcar → median → boxcar → Gaussian → GMT surface 重建。
滤波后的场返回 A 原始坐标和有效覆盖区；边缘外推最多 1.5 个 B 像元间距，超出支持区保持 NaN。

Run5 生成 `iono_correction_filt.grd`、`unwrap_corrected.grd`、校正后包裹相位
`phase_LA_after.grd` 及对应 LOS。Run6 将这个校正后包裹相位送入 SNAPHU，再次解缠并换算最终 LOS，
不重复高斯和 Goldstein。

## 5. 可选的多边形 mask：先预览再计算

普通 run5 不加额外多边形。需要剔除断开区域时，在 run4 后使用掩膜分支。
顶点按多边形边界顺序填写，单位是产品 CRS 的米；`exclude` 剔除圈内，`keep` 只保留圈内。
下列坐标是开发时 T112A 场景的示例，其他数据应重新选点。

```bash
./run5mask_iono_GSLC_LA_LB.csh 1 --mask-mode exclude \
  --mask-polygon 505000 3505000 750000 3575000 \
                 765000 3525000 520000 3455000
```

检查 `iono_correction_mask/preview/LA/mask_preview.pdf` 和对应 LB 预览。
模式 1 保存选区，不计算电离层。确认后执行：

```bash
./run5mask_iono_GSLC_LA_LB.csh 2
./run6_corrected_GSLC_LA.csh \
  --iono-dir iono_correction_mask -o LA_iono_corrected_mask
```

模式 2 读取保存的选区。结果位于 `iono_correction_mask/result/`，经纬度结果在 `result/lonlat/`。
修改选区后先重做模式 1；旧结果移入 `history/`，进程锁阻止同一目录并发运行。
Run6 核对选区与结果，将 mask 延续到校正相位、相干性、解缠和 LOS。

等价的 Python 入口为：

```bash
python gslc_workflow.py mask 1 --mask-mode exclude \
  --mask-polygon 505000 3505000 750000 3575000 \
                 765000 3525000 520000 3455000
# 看过预览后继续。
python gslc_workflow.py mask 2
python gslc_workflow.py correct \
  --iono-dir iono_correction_mask -o LA_iono_corrected_mask
```

## 6. 坐标、LOS 和导出控制

GSLC 已经地理编码。Run4–run6 是以最近邻重投影得到经纬度结果，不走旧 RSLC 的 `trans.dat` 路径。
原生网格保留 CRS；经纬度网格为 EPSG:4326，坐标单位为度。
默认采用 rasterio 内部重投影、netCDF4 写入，不依赖外部 `gdalwarp` 的 HDF5/netCDF 插件。

LOS 按各频段自身波长计算 `-波长/(4*pi) * unwrap`。网格单位为 m，图上为 mm，朝卫星为正。
默认不选参考区、不减常数，所以相位零点和模糊度仍会保留。所有图不加梯度阴影，colorbar 在地图下方。

分步入口支持 `--no-pdf`；run4–run6 此时仍导出经纬度网格，也可加 `--no-lonlat` 关闭经纬度输出。
一键总控支持 `--no-pdf`。单独导出经纬度或可选 GeoTIFF：

```bash
python reproject_gslc.py LA/intf/<pair> --resume --pdf
python reproject_gslc.py LA/intf/<pair> --resume --geotiff
```

日志和恢复见[服务器说明](SERVER_RUN.zh-CN.md)，来源与解释限制见
[UPSTREAM_PROVENANCE.md](UPSTREAM_PROVENANCE.md)。
