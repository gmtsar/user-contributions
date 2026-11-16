# 服务器运行与恢复

[English](SERVER_RUN.md) | **简体中文** · [首页](README.zh-CN.md)

## 安装仓库副本

把 `nisar_gslc/` 包复制到服务器，再进入该包目录。输入数据和处理结果与源码分开，
`.gitignore` 已排除常见数据及生成文件。

激活已有 GMTSAR/SAR 环境；Python 依赖尚未安装时再补齐：

```bash
python -m pip install -r requirements.txt
python -c 'import numpy,scipy,h5py,netCDF4,pyproj,rasterio; print(pyproj.CRS.from_epsg(32646)); print(rasterio.crs.CRS.from_epsg(4326))'
command -v csh phasefilt snaphu gmt gs
export NISAR_GSLC_PYTHON="$(command -v python)"
chmod +x *.csh
```

工具应能在 `PATH` 中找到。分步入口在对应计算步骤支持 `--phasefilt`、`--snaphu`、`--gmt` 指定路径；
一键总控只接受其帮助中列出的选项，不自动转发任意计算参数。

默认经纬度导出由 rasterio 与 netCDF4 内部完成，不依赖外部 `gdalwarp` 的 HDF5/netCDF 驱动。
`pdftoppm` 仅用于可选的 PNG 预览，不影响 GRD/PDF。

## 在仓库目录直接运行

在 `raw/` 放入且仅放入两景 GSLC `.h5` 或软链接。新项目先检查 `configs/` 模板，然后执行：

```bash
nohup ./run1-6_gslc_all.csh > run1-6_gslc_all.log 2>&1 &
tail -f run1-6_gslc_all.log
```

如果想先修改项目配置：

```bash
./run1_prepare_GSLC_LA_LB.csh
# 编辑 LA/config_GSLC_A.txt、LB/config_GSLC_B.txt、config_IONO.txt。
nohup ./run1-6_gslc_all.csh --from-step 2 > run2-6.log 2>&1 &
```

## 将数据项目放在仓库外

将示例路径改成实际位置，在该项目 `raw/` 下放两景数据。以下命令从代码目录运行：

```bash
export GSLC_PROJECT=/path/to/data_project
mkdir -p "$GSLC_PROJECT/raw"
# 先把两景 GSLC 或软链接放入 "$GSLC_PROJECT/raw"。
./run1_prepare_GSLC_LA_LB.csh --project "$GSLC_PROJECT"
# 编辑 "$GSLC_PROJECT" 下生成的配置。
nohup ./run1-6_gslc_all.csh --project "$GSLC_PROJECT" \
  --from-step 2 > run2-6.log 2>&1 &
```

各个分步入口或直接 Python 调用也必须使用相同的 `--project`，例如：

```bash
python gslc_workflow.py process --project "$GSLC_PROJECT"
```

同一个项目的分步计算和一键总控不要同时启动。

## 日志与最终结果

外层日志显示每步开始、完成和失败，并打印详细日志目录
`PROJECT/logs/all_UTC时间戳_PID/`。具体计算查看其中的 `run2.log`、`run4.log` 等。
某一步失败后，总控停止。

默认最终结果为：

```text
PROJECT/LA_iono_corrected/unwrap.grd
PROJECT/LA_iono_corrected/los.grd
PROJECT/LA_iono_corrected/unwrap.pdf
PROJECT/LA_iono_corrected/los.pdf
PROJECT/LA_iono_corrected/lonlat/unwrap_ll.grd
PROJECT/LA_iono_corrected/lonlat/los_ll.grd
```

A/B 顺序处理，单频段内按 `unwrap_nproc` 并行求解 SNAPHU 块。
整幅运行时间取决于覆盖范围、相干性、网格间距与服务器资源，不给固定时长承诺。

## 失败后继续

总控支持 `--from-step N`（N=1…6），从指定步骤一直运行到第六步。
已有输出仍受保护，该选项不是通用覆盖开关。

原生解缠已完成、经纬度导出失败时，可以：

```bash
nohup ./run1-6_gslc_all.csh --from-step 4 --resume \
  > continue_run4-6.log 2>&1 &
```

`--resume` 作用于 run4，核对完成报告和参数，复用完整的解缠/LOS，保留已记录失败的经纬度目录后重试，
再继续未完成的频段及 run5/run6。它不恢复未完成的 SNAPHU 子块，也不自动覆盖已有 run5/run6 结果。

只恢复 run4 或只重试经纬度导出：

```bash
./run4_unwrap_GSLC_LA_LB.csh --resume
python reproject_gslc.py LA/intf/<pair> --resume --pdf
```

将 `<pair>` 替换为实际目录。外部数据项目按需使用 `--project` 或对应的完整数据路径。

想对比另一组电离层滤波参数时，修改 `config_IONO.txt`，并使用新输出目录：

```bash
./run5_iono_GSLC_LA_LB.csh -o iono_correction_trial
./run6_corrected_GSLC_LA.csh \
  --iono-dir iono_correction_trial -o LA_iono_corrected_trial
```

可选多边形分支使用已保存的预览和模式 2 结果，见[流程说明](WORKFLOW_GSLC.zh-CN.md)。
直接调用 Python 的方式见 [README](README.zh-CN.md)。

## 验证范围

复制的源码对应 v0.6.5。本地小范围诊断已检查 A/B 处理、掩膜校正、经纬度导出恢复及重复图复用。
A/B 有效视数仍为临时设置。新场景的整幅耗时和科学有效性，应按自身结果评估；
解释限制见 [README](README.zh-CN.md)。
