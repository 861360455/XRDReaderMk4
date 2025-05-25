import streamlit as st
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
import itertools
import re
from itertools import combinations
# 引入 pymatgen 的 XRDCalculator
from mp_api.client import MPRester
from pymatgen.analysis.diffraction.xrd import XRDCalculator

# -------------------- 新版 fetch_materials --------------------
def get_chemsys(elements):
    """
    根据元素列表返回化学系统字符串（例如 "C-Ni"）。
    """
    return "-".join(sorted(elements))

def fetch_materials(allowed_elements, api_key):
    """
    根据允许的元素生成所有非空子集，并利用 Materials Project 查询对应材料。
    返回每个记录包含 material_id、formula、structure、cif 字符串，
    同时附加 label 字段（格式为 "material_id (formula)"）。
    """
    all_subsets = []
    for r in range(1, len(allowed_elements) + 1):
        all_subsets.extend(itertools.combinations(allowed_elements, r))
    
    results = []
    with MPRester(api_key) as mpr:
        for subset in all_subsets:
            chemsys_str = get_chemsys(subset)
            st.write(f"查询化学系统: {chemsys_str}")
            docs = mpr.materials.search(
                chemsys=chemsys_str, 
                fields=["structure", "material_id", "formula_pretty"]
            )
            for doc in docs:
                cif_str = doc.structure.to(fmt="cif")
                label = f"{doc.material_id} ({doc.formula_pretty})"
                results.append({
                    "chemsys": chemsys_str,
                    "material_id": doc.material_id,
                    "formula": doc.formula_pretty,
                    "structure": doc.structure,
                    "cif": cif_str,
                    "label": label
                })
                st.write(f"  Material ID: {doc.material_id}, Formula: {doc.formula_pretty}")
    return results

def step1_fetch_materials():
    st.subheader("Step 1: 输入材料体系 & 获取 MP 数据")
    if "step1_done" not in st.session_state:
        st.session_state["step1_done"] = False
    if "materials" not in st.session_state:
        st.session_state["materials"] = None

    elements_input = st.text_input("化学系统元素（空格分隔）", value="C Ni")
    api_key_input = st.text_input("Materials Project API key", value="vfxI2hPnzLupiP4cdjuSQDBlEg3KjC1Q")

    if st.button("Run Step 1"):
        allowed_elements = elements_input.split()
        api_key = api_key_input.strip()
        st.write("开始从 Materials Project 查询...")
        materials = fetch_materials(allowed_elements, api_key)  # 假设该函数已定义
        st.session_state["materials"] = materials
        st.write(f"共获得 {len(materials)} 个材料记录。")
        # 从材料数据中提取 label 列表供后续使用
        material_labels = [m["label"] for m in materials]
        st.session_state["material_labels"] = material_labels
        st.session_state["step1_done"] = True
        # Step1 结束后
        material_formulas = [m["formula"] for m in materials]
        st.session_state["material_formulas"] = material_formulas
        st.success("Step 1 完成：已获取 MP 数据")

# -------------------- Step 2: 计算理论离散峰 --------------------
def step2_compute_theoretical_xrd():
    st.subheader("Step 2: 设置发射源 & 计算理论离散峰")
    if not st.session_state.get("step1_done", False):
        st.warning("请先完成 Step 1。")
        return
    if "step2_done" not in st.session_state:
        st.session_state["step2_done"] = False

    # 选择 XRD 发射源
    st.write("选择或指定 XRD 发射源：")
    xrd_source_options = ["CuKa", "MoKa", "CoKa", "CrKa", "自定义(输入波长)"]
    chosen_source = st.selectbox("XRD 发射源", xrd_source_options, index=0)
    if chosen_source == "CuKa":
        wavelength = 1.5406
    elif chosen_source == "MoKa":
        wavelength = 0.7093
    elif chosen_source == "CoKa":
        wavelength = 1.78897
    elif chosen_source == "CrKa":
        wavelength = 2.2897
    else:
        wavelength = st.number_input("自定义波长 (Å)", value=1.5406, step=0.01)

    # -------------------- Step 2: 计算理论离散峰 --------------------
    if st.button("Run Step 2: 计算理论离散峰"):
        st.write("开始计算理论离散峰...")
        materials = st.session_state["materials"]  # 假设已经在上一步存好的
        discrete_patterns = []
        material_formulas = st.session_state["material_formulas"]  # 由 Step1 存的化学式列表
        materials = st.session_state["materials"]
        # 创建 XRDCalculator
        wavelength = 1.5418  # 例如 Cu Kα (Angstrom)，用户可根据需求修改
        calc = XRDCalculator(wavelength=wavelength)
        
        # 对每个材料结构直接计算离散衍射峰
        for mat in materials:
            pattern = calc.get_pattern(mat["structure"])
            discrete_patterns.append(pattern)

        # 存到 session
        st.session_state["discrete_patterns"] = discrete_patterns
        st.session_state["step2_done"] = True
        st.success("Step 2 完成：理论离散峰计算完成。")


def step3_peak_matching():
    st.subheader("Step 3: 多峰匹配及标注")
    if not st.session_state.get("step2_done", False):
        st.warning("请先完成 Step 2。")
        return
    
    # 如果 session 中尚未有 locked_phases 字典，就初始化一下
    if "locked_phases" not in st.session_state:
        st.session_state["locked_phases"] = {}
    
    # 1. 读取实验谱
    uploaded_file = st.file_uploader("上传实验 XRD 谱文件", type=["xy", "csv"])
    if not uploaded_file:
        st.info("请先上传实验 XRD 文件")
        return

    # 读取实验数据
    try:
        # 若文件是空格分隔或逗号分隔，需要酌情修改下面的参数
        df = pd.read_csv(uploaded_file, delim_whitespace=True, header=None, names=["two_theta", "intensity"])
        st.write("实验谱预览：", df.head())
        exp_x = df["two_theta"].values
        exp_y = df["intensity"].values
        st.session_state["exp_x"] = exp_x
        st.session_state["exp_y"] = exp_y
    except Exception as e:
        st.error(f"读取 XRD 文件出错：{e}")
        return

    # 2. 获取 Step2 生成的 pattern
    discrete_patterns = st.session_state["discrete_patterns"]
    material_formulas = st.session_state["material_formulas"]  # 存在 session_state 中
    def pick_top3_peaks_with_merge(peaks, tol=0.2):
        """
        给定离散XRD峰列表( (two_theta, intensity) )，先取强度最大的3条峰，
        如果它们之间有两条相距 < 2*tol，则合并该对峰并重新排序、再取3条峰重复检测。
        直至不再需要合并，返回最终的3条峰（若实际峰数<3则全部返回）。
        
        合并规则：
        - 新峰坐标 = 两峰中强度更大那条峰的坐标
        - 新峰强度 = 两峰中更大的强度
        """
        if not peaks:
            return []

        # 拷贝一份，避免对原列表修改
        current_peaks = list(peaks)

        while True:
            # 按强度降序排序
            current_peaks.sort(key=lambda x: x[1], reverse=True)
            
            # 如果峰不足3个，就直接返回全部
            if len(current_peaks) < 3:
                return current_peaks
            
            # 先截取当前最强的3个
            top3 = current_peaks[:3]

            # 检查 top3 中是否有两条峰的距离 < 2*tol
            merged = False
            for peakA, peakB in combinations(top3, 2):
                if abs(peakA[0] - peakB[0]) < 2*tol:
                    # 准备合并
                    if peakA[1] >= peakB[1]:
                        new_peak = (peakA[0], peakA[1])
                    else:
                        new_peak = (peakB[0], peakB[1])

                    # 在 current_peaks 中移除这两个峰，再加入合并后的新峰
                    current_peaks.remove(peakA)
                    current_peaks.remove(peakB)
                    current_peaks.append(new_peak)
                    
                    # 标记这一次合并完成，结束 for 循环并进入下一轮 while
                    merged = True
                    break
            
            # 如果这一轮没有发生任何合并，就说明已满足条件，可以退出循环
            if not merged:
                break

        # 最后再按强度降序排一下，取最强3条
        current_peaks.sort(key=lambda x: x[1], reverse=True)
        return current_peaks[:3]
    # 构建三强峰字典
    new_candidate_dict = {}
    for formula, pattern in zip(material_formulas, discrete_patterns):
        all_peaks = list(zip(pattern.x, pattern.y))
        
        # 调用我们刚刚定义的选取+合并函数
        top3_peaks = pick_top3_peaks_with_merge(all_peaks, tol=0.2)
        
        new_candidate_dict[formula] = top3_peaks

    # 3. 实验峰检测
    exp_peaks_idx, exp_props = find_peaks(exp_y, prominence=10, distance=20)
    exp_peaks_all = exp_x[exp_peaks_idx]
    # 按 prominence 降序排序
    sorted_indices = np.argsort(exp_props["prominences"])[::-1]
    exp_peaks_sorted = exp_peaks_all[sorted_indices]

    if exp_peaks_sorted.size == 0:
        st.error("未检测到实验谱峰，请检查实验数据！")
        return

    n_peak = st.slider("可视化时显示多少个峰（从最显著开始）",
                       min_value=1, max_value=int(len(exp_peaks_sorted)),
                       value=int(len(exp_peaks_sorted)), step=1)
    exp_peaks_to_show = exp_peaks_sorted[:n_peak]

    # 系统性偏移（滑条）
    offset_min = -5.0
    offset_max = 5.0
    offset_step = 0.1
    offset = st.slider("系统性偏移修正 (°)", min_value=offset_min, max_value=offset_max, value=0.0, step=offset_step)
    
    # 找到要显示的峰在 exp_peaks_all 中的索引（这些是排序后的）
    peaks_to_show_idx_in_all = sorted_indices[:n_peak]
    # 获取这些峰对应的原始索引（在 exp_x/exp_y 中的索引）
    exp_peaks_idx_to_show = exp_peaks_idx[peaks_to_show_idx_in_all]

    # 构建只包含可视化显示的峰的列表（偏移前和偏移后）
    exp_peaks_data = [(exp_x[i], exp_y[i]) for i in exp_peaks_idx_to_show]
    exp_peaks_data_corrected = [(p[0] + offset, p[1]) for p in exp_peaks_data]

    # ④ 多峰匹配逻辑
    tol = 0.2       # 2θ 容差
    ratio_tol = 0.3 # 相对强度比 容差 (±30% 仅供演示)
    
    identified_phases = []  
    matched_peaks_dict = {}  # 存储每个相匹配成功的实验峰下标，用于后面标注


    # (2) 在匹配阶段，不再跳过 <3 条峰的相
    #     原本是： if len(top3_theory_peaks) < 3: continue
    #     现在改为：任何数量都参与匹配
    identified_phases = []
    matched_peaks_dict = {}

    for phase, theory_peaks in new_candidate_dict.items():
        peak_count = len(theory_peaks)
        if peak_count == 0:
            # 没有峰就确实无法匹配，跳过
            continue

        # =============== 逐条理论峰做匹配 ===============
        match_count = 0
        matched_indices_for_this_phase = []

        for (tx, ty) in theory_peaks:
            diffs = [abs(e[0] - tx) for e in exp_peaks_data_corrected]
            min_idx = np.argmin(diffs)
            if diffs[min_idx] <= tol:
                match_count += 1
                matched_indices_for_this_phase.append(min_idx)
            else:
                # 只要有一条峰匹配失败，则可以判定为不完全匹配
                break

        # =============== 判断是否成功匹配 ===============
        # 如果全部峰都匹配成功，则进一步判断强度比（若有需要）
        if match_count == peak_count:
            # 只有当峰数 >=2 时，才有必要做强度比
            if peak_count >= 2:
                theo_ints = [p[1] for p in theory_peaks]
                theo_max_int = max(theo_ints)
                theo_min_int = min(theo_ints)
                # 避免除0
                if theo_min_int == 0:
                    continue

                exp_ints_matched = [exp_peaks_data_corrected[i][1]
                                    for i in matched_indices_for_this_phase]
                if min(exp_ints_matched) == 0:
                    continue

                theo_ratio = theo_max_int / theo_min_int
                exp_ratio = max(exp_ints_matched) / min(exp_ints_matched)
                ratio_diff = abs(exp_ratio - theo_ratio) / theo_ratio

                if ratio_diff <= ratio_tol:
                    identified_phases.append(phase)
                    matched_peaks_dict[phase] = matched_indices_for_this_phase
            else:
                # 只有1条峰可匹配的情况，直接视为通过
                identified_phases.append(phase)
                matched_peaks_dict[phase] = matched_indices_for_this_phase

    # ⑤ 可视化
    fig, ax = plt.subplots(figsize=(10, 6))
    # 先画整体实验谱
    ax.plot(exp_x, exp_y, label="实验谱")

    # 用红叉标示"n_peak"个最显著峰 (只是作参考，可改成全部标示)
    exp_peaks_corrected_to_show = exp_peaks_to_show + offset
    for peak_val in exp_peaks_to_show:
        idx = np.argmin(np.abs(exp_x - peak_val))
        ax.plot(exp_x[idx], exp_y[idx], "rx")

    # 改进版：在每个峰上动态垂直堆叠标注，避免重叠
    peak_annotation_count = dict()  # key: corrected x, value: number of annotations stacked

    # 1️⃣ 先绘制锁定的相（不受滑条影响，始终显示）
    for phase, lock_info in st.session_state["locked_phases"].items():
        # 取出已经保存好的 peak 坐标
        peaks_to_plot = lock_info["peaks"]
        peak_color = "green"
        label_prefix = "🔒 "
        for px, py in peaks_to_plot:
            key = round(px, 4)
            count = peak_annotation_count.get(key, 0)
            ax.plot(px, py, "o", color=peak_color)
            ax.annotate(f"{label_prefix}{phase}",
                        (px, py),
                        textcoords="offset points",
                        xytext=(0, 8 + 10 * count),
                        ha="center",
                        color=peak_color,
                        fontsize=8)
            peak_annotation_count[key] = count + 1

    # 2️⃣ 再绘制当前偏移量下动态匹配的相
    for phase, idx_list in matched_peaks_dict.items():
        if phase in st.session_state["locked_phases"]:
            continue  # 锁定的相已经画过了，避免重复
        matched_indices = idx_list
        peaks_to_plot = [exp_peaks_data[i] for i in matched_indices]
        peak_color = "blue"
        label_prefix = ""
        for px, py in peaks_to_plot:
            key = round(px, 4)
            count = peak_annotation_count.get(key, 0)
            ax.plot(px, py, "o", color=peak_color)
            ax.annotate(f"{label_prefix}{phase}",
                        (px, py),
                        textcoords="offset points",
                        xytext=(0, 8 + 10 * count),
                        ha="center",
                        color=peak_color,
                        fontsize=8)
            peak_annotation_count[key] = count + 1

    ax.set_xlabel("2θ (°)")
    ax.set_ylabel("Intensity")
    ax.set_title("实验谱（多峰匹配标注）")
    ax.legend()
    st.pyplot(fig)

    # 显示匹配结果，并支持锁定/解锁
    st.write("匹配到的相：")
    if len(identified_phases) == 0:
        st.warning("未找到符合三强峰匹配条件的相。")
    else:
        for phase in identified_phases:
            col1, col2 = st.columns([3, 1])
            with col1:
                st.write(f"- {phase} (三强峰均匹配，强度比通过)")
            with col2:
                if phase not in st.session_state["locked_phases"]:
                    if st.button(f"🔒 锁定 {phase}", key=f"lock_{phase}"):
                        matched_indices_for_this_phase = matched_peaks_dict[phase]
                        peaks_to_lock = [exp_peaks_data[i] for i in matched_indices_for_this_phase]
                        st.session_state["locked_phases"][phase] = {
                            "peaks": peaks_to_lock,
                            "offset": offset,
                            "matched_indices": matched_indices_for_this_phase
                        }
                else:
                    if st.button(f"🔓 解锁 {phase}", key=f"unlock_{phase}"):
                        del st.session_state["locked_phases"][phase]
                    else:
                        st.markdown("✅ 已锁定")

    # 显示清除所有锁定按钮
    if st.session_state["locked_phases"]:
        if st.button("🗑️ 清除所有锁定的相"):
            st.session_state["locked_phases"].clear()
            st.success("已清除所有锁定的相。")

    # 用于生成给定区间内的对称 offset 序列（包含 0.0）
    def generate_offset_sequence(start, stop, step):
        offsets = [0.0]
        val = step
        while val <= stop:
            offsets.append(-val)
            offsets.append(val)
            val += step
        # 排序并去重，最后返回区间内有效的 offset
        offsets = sorted(set(offsets))
        return [o for o in offsets if start <= o <= stop]

    # 自动批量匹配
    if st.button("🚀 自动批量匹配并锁定所有相"):
        already_locked = set(st.session_state["locked_phases"].keys())
        offset_list = generate_offset_sequence(offset_min, offset_max, offset_step)
        auto_matched_count = 0

        for offset_try in offset_list:
            # 逐个相尝试匹配
            for phase, top3_theory_peaks in new_candidate_dict.items():
                if phase in already_locked:
                    continue
                if len(top3_theory_peaks) < 3:
                    continue

                # 理论强度比
                theo_ints = [p[1] for p in top3_theory_peaks]
                if min(theo_ints) == 0:
                    continue
                theo_ratio = max(theo_ints) / min(theo_ints)

                # 当前 offset 下的实验峰坐标
                exp_peaks_data_corrected_try = [(p[0] + offset_try, p[1]) for p in exp_peaks_data]

                match_count = 0
                matched_indices = []
                # 尝试匹配三强峰
                for tx, ty in top3_theory_peaks:
                    diffs = [abs(e[0] - tx) for e in exp_peaks_data_corrected_try]
                    min_idx = np.argmin(diffs)
                    if diffs[min_idx] <= tol:
                        match_count += 1
                        matched_indices.append(min_idx)
                    else:
                        break

                if match_count == 3:
                    matched_exp_ints = [exp_peaks_data_corrected_try[idx][1] for idx in matched_indices]
                    if min(matched_exp_ints) == 0:
                        continue
                    exp_ratio = max(matched_exp_ints) / min(matched_exp_ints)
                    ratio_diff = abs(exp_ratio - theo_ratio) / theo_ratio
                    if ratio_diff <= ratio_tol:
                        # 匹配成功，自动锁定
                        peaks_to_lock = [(exp_peaks_data[idx][0], exp_peaks_data[idx][1]) 
                                         for idx in matched_indices]
                        st.session_state["locked_phases"][phase] = {
                            "peaks": peaks_to_lock,
                            "offset": offset_try,
                            "matched_indices": matched_indices
                        }
                        auto_matched_count += 1
                        already_locked.add(phase)

        st.success(f"自动匹配并锁定完成，共锁定 {auto_matched_count} 个相。")
    st.markdown("### 🔒 已锁定相")
    locked_phases = st.session_state.get("locked_phases", {})
    if not locked_phases:
        st.info("当前没有已锁定的相。")
    else:
        for phase, info in locked_phases.items():
            offset_val = info["offset"]
            st.write(f"- `{phase}` 锁定于偏移量 `{offset_val:.2f}°`")
    if st.button("🔍 显示锁定相的三强峰原始位置"):
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.plot(exp_x, exp_y, label="实验谱")

        for phase, info in st.session_state["locked_phases"].items():
            offset_used = info["offset"]
            indices = info["matched_indices"]
            
            # 匹配到的实验峰位置（已偏移）
            matched_exp_positions = [exp_x[i] + offset_used for i in indices]
            
            # 找这个相的原始三强理论峰
            if phase in new_candidate_dict:
                top3_theory_peaks = new_candidate_dict[phase]
                
                # 反向修正坐标：减去 offset，还原三强峰原始坐标
                corrected_theory_positions = [x - offset_used for x, _ in top3_theory_peaks]
                
                # 绘制三强峰原始坐标
                for px in corrected_theory_positions:
                    ax.axvline(px, color='orange', linestyle='--', linewidth=1)
                    ax.annotate(f"{phase}", (px, max(exp_y) * 0.7),
                                rotation=90, color="orange", fontsize=8,
                                xytext=(0, 5), textcoords="offset points", ha="center")
                # 2️⃣ 同时在界面打印它们的 2θ 坐标
                st.write(f"**锁定相：{phase}**")
                st.write("三强峰原始(未偏移) 2θ 坐标：")
                for p in corrected_theory_positions:
                    st.write(f" - {p:.3f}°")
                st.write("---")  # 分割线，方便查看
        ax.set_title("锁定相的三强峰（原始未偏移坐标）")
        ax.set_xlabel("2θ (°)")
        ax.set_ylabel("Intensity")
        st.pyplot(fig)
# -------------------- 主程序入口 --------------------
def main():
    st.title("自动化 XRD 数据处理与峰匹配")
    step = st.sidebar.selectbox("选择步骤", ["Step 1: 获取材料", "Step 2: 计算理论谱", "Step 3: 峰匹配"])
    if step == "Step 1: 获取材料":
        step1_fetch_materials()
    elif step == "Step 2: 计算理论谱":
        step2_compute_theoretical_xrd()
    elif step == "Step 3: 峰匹配":
        step3_peak_matching()

if __name__ == "__main__":
    main()
