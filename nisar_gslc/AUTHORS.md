# Authors and Acknowledgments

[English](#english) | [简体中文](#中文说明)

## English

### GSLC workflow author and maintainer

- **Name:** Xin Wang
- **Affiliation:** University of Science and Technology of China (USTC)
- **Location:** Hefei, China
- **Contact:** [xinw11@mail.ustc.edu.cn](mailto:xinw11@mail.ustc.edu.cn)
- **GitHub:** [Xin-Wang-520](https://github.com/Xin-Wang-520)
- **Project:** [gmtsar-nisar-gslc-dinsar-workflow](https://github.com/Xin-Wang-520/gmtsar-nisar-gslc-dinsar-workflow)

This attribution covers the GSLC workflow adaptation, integration and
maintenance: map-grid processing, configuration and staged entry points,
A/B workflow integration, geographic exports, plotting and documentation.

### Upstream methods and code

The ionospheric frequency combination follows Gomba et al. (2016), as
implemented in the earlier
[RSLC workflow](https://github.com/Xin-Wang-520/gmtsar-nisar-dinsar-workflow).
That workflow credits subsequent script development to **Yajun Zhang** and
NISAR LA/LB workflow integration to **Xin Wang**.

`nearest_grid_c_exact.py` is reproduced from the recorded upstream revision;
its existing Xin Wang author header is retained. External software, including
GMTSAR, GMT and SNAPHU, retains its respective authorship and licenses.
Package authorship does not replace these method and software credits.

See [UPSTREAM_PROVENANCE.md](UPSTREAM_PROVENANCE.md) for the source revision,
adaptation details and attribution. The existing distribution license is
provided in [LICENSE](LICENSE).

## 中文说明

### GSLC 工作流作者与维护者

- **姓名：** Xin Wang
- **单位：** 中国科学技术大学（USTC）
- **所在地：** 中国合肥
- **联系邮箱：** [xinw11@mail.ustc.edu.cn](mailto:xinw11@mail.ustc.edu.cn)
- **GitHub：** [Xin-Wang-520](https://github.com/Xin-Wang-520)
- **项目：** [gmtsar-nisar-gslc-dinsar-workflow](https://github.com/Xin-Wang-520/gmtsar-nisar-gslc-dinsar-workflow)

上述署名对应 GSLC 工作流的适配、集成与维护，包括地图网格处理、配置与分步入口、
A/B 流程集成、经纬度导出、绘图和文档。

### 上游方法与代码致谢

电离层频率组合沿用原 [RSLC 工作流](https://github.com/Xin-Wang-520/gmtsar-nisar-dinsar-workflow)
实现的 Gomba 等（2016）方法。该工作流中，后续脚本开发署名为 **Yajun Zhang**，
NISAR LA/LB 流程集成署名为 **Xin Wang**。

`nearest_grid_c_exact.py` 来自已记录的上游版本，保留原有 Xin Wang 作者头。
GMTSAR、GMT、SNAPHU 等外部软件保留各自作者与许可证信息。
本项目署名不替代上述方法和软件的来源署名。

来源版本、适配细节与署名见 [UPSTREAM_PROVENANCE.md](UPSTREAM_PROVENANCE.md)，
现有发布许可证见 [LICENSE](LICENSE)。
