# ADLFC 电力系统负荷频率控制仿真

对应论文 *A Novel Adaptive Data-Driven Load Frequency Control for Unknown Nonlinear Power Systems* 的数值仿真代码。

本包采用已与当前论文图片逐文件核准的版本：A 组为 2026-10-23 取消控制输入幅值截断后的三种噪声实验；B 组为 2026-10-27 在 60、200、360 s 分时接入阶跃、多正弦、OU 随机负荷的三方法对比实验。核心模型和控制器文件原样保留。

## 运行

已在 Windows、MATLAB R2023b 验证。将 MATLAB 当前目录切换到仓库根目录，执行：

```matlab
run_all
```

依次运行 A/B 六个算例、导出 Tables III–V、进行参考数据核验，并输出 11 张仿真 PNG 图。不需要原论文项目中的历史文件夹，不需要 Simulink 或 Python。

也可单独执行 `run_experiment_a`、`run_experiment_b`、`compute_paper_metrics`、`verify_results`、`plot_results`。离线采集与权重生成可单独运行 `[w, certificate] = compute_adlfc_weights;`。

## 重要设置

- A/B 仿真时长分别为 600/420 s，均采用 0.001 s 前向欧拉积分。
- 各区域滤波时间常数实际为 0.01、0.01、0.02 s。
- 调速器完整导数限制为 ±2 p.u./s；控制输入不施加额外幅值截断。
- B 组三方法每个积分时刻采用相同负荷。OU 种子为 20261027，多正弦和 OU 采用 5 s 余弦平滑接入。
- 默认读取本次作图实际使用的固定权重；离线采样矩阵、设计中间量和诊断记录一并提供。
- 频差直接使用 rad/s 计算 ACE，详细公式、偏置、基线增益及初值见英文 README、配置文件和源代码。

## 输出与核验

所有新生成的全分辨率结果、CSV 指标表、核验报告与图片均写入 `results/`，该目录不提交到 GitHub。`verify_results` 对照原始六次作图运行的参考采样及稿件精度下的三张指标表。

图片使用本次重算轨迹，保留抽样块中的真实极值；仓库绘图采用简洁布局，与论文排版版式不同。

完整参数和初值保存在 `config/frozen_settings.mat` 及可读的 JSON 中。来源文件与 SHA256 见 `docs/source_provenance.json`。

按作者要求，当前暂不添加许可证。
