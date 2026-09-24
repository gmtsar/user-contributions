# xcorr3 详细用法

[返回首页](README.zh-CN.md)

将 PATH 设置写入 shell 启动文件即可保留。MT 默认使用 Ubuntu x86-64 库路径，
其他环境需调整后端 Makefile。CUDA 可设置 `CUDA_HOME`、`CUDA_ARCH`，
例如 `make cc CUDA_ARCH=sm_75`。CC 的 `-geocode` 额外需要 GMT、GMTSAR
`proj_ra2ll.csh`、C shell 及匹配的 `trans.dat`，默认不启用。

## 使用示例

```sh
xcorr3 master.PRM secondary.PRM [参数]
xcorr3 master.PRM secondary.PRM --backend mt -nproc 6 -nx 20 -ny 50
xcorr3 --backend cc master.PRM secondary.PRM -xsearch 128 -ysearch 128
xcorr3 --help
```

不足两份输入文件时显示所需参数、用法和示例，并返回退出码 2。
两份 PRM 记录的 SLC 路径应能从当前工作目录访问。

| 参数 | 用途 |
|---|---|
| `--backend original\|mt\|cc` | 选择后端，不自动降级切换 |
| `-nproc N` / `--nproc N` | MT 进程数，或 `auto`；命令行优先 |
| `-nx N -ny N` | 距离向、方位向采样点数 |
| `-freq` | 频域互相关；CC 自动映射为其默认模式 |
| `-xsearch N -ysearch N` | 搜索半径 |
| `-range_interp N -interp N` | 插值倍数 |
| `-noshift -norange -nointerp` | 忽略粗偏移、禁用相应插值 |
| `--config FILE` | 从现有 GMTSAR 配置读取后端及进程数 |
| `--dry-run` | 仅显示选择及参数，不处理数据 |

在已有配置中各追加一次 `xcorr_backend = mt` 和 `xcorr_nproc = 6`，
参见 [配置示例](xcorr-options.config)。命令行覆盖配置；缺省为原版和自动进程数。
配置中的进程数仅对 MT 生效。GMTSAR 自身不会读取这些新增选项。

```sh
xcorr3 --config config.txt master.PRM secondary.PRM -nx 20 -ny 50
# 可选：转接前台、非 TOPS 流程中的裸 xcorr 调用。
xcorr3 --config config.txt --run p2p_processing.csh ALOS a b config.txt
```

`--run` 只为子进程设置临时 PATH；绝对路径调用或重置 PATH 的脚本会绕过转接。
原版 TOPS 几何配准/ESD 不调用 xcorr，不能由本入口加速替代。
被转接的 xcorr 失败后，即使外层脚本忽略错误并返回成功，入口仍返回非零。

保留各后端输出格式：通常为 `freq_xcorr.dat`，原版/MT 五列，CC 六列（新增
`peak_snr`）。CC 不支持时域和实数/网格输入模式，结果可能与 CPU 不同。
同一输出目录只运行一个任务，使用输出前须检查退出码。

[测试与验证](VALIDATION.md) · [MT 说明](backends/mt/README.md) ·
[CC 说明](backends/cc/README.zh-CN.md)。PowerShell 可运行 Python 帮助和入口测试；
科学计算后端及 `--run` 仍要求 Linux/POSIX。
