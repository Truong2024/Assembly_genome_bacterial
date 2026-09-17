#!/bin/bash
# chmod +x myscript_ver2.sh && ./myscript_ver2.sh          -> hiện menu để chọn
# hoặc: ./myscript_ver2.sh 4                                -> chạy thẳng bước 4 (không cần menu)
# Pipeline dạng menu, có checkpoint: bỏ qua bước đã hoàn thành.

set -e

# ══════════════════════════════════════════════════════════
# MÀU SẮC & HÀM LOG
# ══════════════════════════════════════════════════════════
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
CYAN='\033[0;36m'
RED='\033[0;31m'
BOLD='\033[1m'
NC='\033[0m'

log_skip()  { echo -e "${YELLOW}[SKIP]${NC}  $1 – đã tồn tại, bỏ qua."; }
log_run()   { echo -e "${GREEN}[RUN]${NC}   $1"; }
log_info()  { echo -e "${CYAN}[INFO]${NC}  $1"; }
log_err()   { echo -e "${RED}[LỖI]${NC}  $1"; }
log_done()  { echo -e "${GREEN}${BOLD}[XONG]${NC} $1"; }

# ══════════════════════════════════════════════════════════
# HÀM TIỆN ÍCH
# ══════════════════════════════════════════════════════════
exists() {
    local target="$1"
    if [ -d "$target" ]; then
        [ "$(ls -A "$target" 2>/dev/null)" ] && return 0 || return 1
    fi
    [ -s "$target" ] && return 0 || return 1
}

env_exists() {
    conda env list | awk '{print $1}' | grep -qx "$1"
}

tool_in_env() {
    conda run -n "$1" which "$2" &>/dev/null
}

init_dirs() {
    mkdir -p QC_report assembly_fastq assembly busco_result quast_result busco_download
}

# ══════════════════════════════════════════════════════════
# -1: CÀI ĐẶT MÔI TRƯỜNG CONDA + CÔNG CỤ
# ══════════════════════════════════════════════════════════
step_install() {
    if ! command -v conda &>/dev/null; then
        log_err "Không tìm thấy 'conda'. Cài Miniconda/Anaconda và thêm vào PATH trước."
        exit 1
    fi

    # pipeline_env: fastqc, fastp, quast
    if env_exists "pipeline_env"; then
        log_skip "Môi trường 'pipeline_env'"
    else
        log_run "Tạo môi trường 'pipeline_env' (fastqc, fastp, quast)"
        conda create -y -n pipeline_env -c bioconda -c conda-forge fastqc fastp quast
    fi
    for tool_spec in "fastqc:fastqc" "fastp:fastp" "quast:quast"; do
        bin_name="${tool_spec%%:*}"
        pkg_name="${tool_spec##*:}"
        if tool_in_env "pipeline_env" "$bin_name"; then
            log_skip "Công cụ '$bin_name' trong pipeline_env"
        else
            log_run "Cài '$pkg_name' vào pipeline_env"
            conda install -y -n pipeline_env -c bioconda -c conda-forge "$pkg_name"
        fi
    done

    # unicycler_env
    if env_exists "unicycler_env"; then
        log_skip "Môi trường 'unicycler_env'"
    else
        log_run "Tạo môi trường 'unicycler_env' (unicycler)"
        conda create -y -n unicycler_env -c bioconda -c conda-forge unicycler
    fi
    if tool_in_env "unicycler_env" "unicycler"; then
        log_skip "Công cụ 'unicycler' trong unicycler_env"
    else
        log_run "Cài 'unicycler' vào unicycler_env"
        conda install -y -n unicycler_env -c bioconda -c conda-forge unicycler
    fi

    # busco_env
    if env_exists "busco_env"; then
        log_skip "Môi trường 'busco_env'"
    else
        log_run "Tạo môi trường 'busco_env' (busco)"
        conda create -y -n busco_env -c bioconda -c conda-forge busco
    fi
    if tool_in_env "busco_env" "busco"; then
        log_skip "Công cụ 'busco' trong busco_env"
    else
        log_run "Cài 'busco' vào busco_env"
        conda install -y -n busco_env -c bioconda -c conda-forge busco
    fi

    log_done "Cài đặt môi trường & công cụ hoàn tất."
}

# ══════════════════════════════════════════════════════════
# 0: ĐỔI TÊN FILE FASTQ.GZ
# ══════════════════════════════════════════════════════════
step_rename() {
    if exists "SRread_1.fastq.gz" && exists "SRread_2.fastq.gz"; then
        log_skip "Đổi tên file (SRread_1/2.fastq.gz đã có)"
    else
        log_run "Đổi tên *R1/R2.fastq.gz -> SRread_1/2.fastq.gz"
        mv *R1.fastq.gz SRread_1.fastq.gz
        mv *R2.fastq.gz SRread_2.fastq.gz
    fi
    log_done "Đổi tên file xong. File: SRread_1.fastq.gz, SRread_2.fastq.gz (thư mục hiện tại)"
}

# ══════════════════════════════════════════════════════════
# 1: FASTQC CHO READ GỐC
# ══════════════════════════════════════════════════════════
step_qc1() {
    if exists "QC_report/SRread_1_fastqc.html" && exists "QC_report/SRread_2_fastqc.html"; then
        log_skip "FastQC read gốc (report đã có)"
    else
        log_run "FastQC cho read gốc"
        conda run --no-capture-output -n pipeline_env fastqc -o QC_report/ SRread_1.fastq.gz
        conda run --no-capture-output -n pipeline_env fastqc -o QC_report/ SRread_2.fastq.gz
    fi
    log_done "FastQC read gốc xong. File: QC_report/SRread_1_fastqc.html, QC_report/SRread_2_fastqc.html"
}

# ══════════════════════════════════════════════════════════
# 2: TRIMMING VỚI FASTP
# ══════════════════════════════════════════════════════════
step_trim() {
    if exists "SRread_1_trimmed.fastq.gz" && exists "SRread_2_trimmed.fastq.gz"; then
        log_skip "Trimming (file trimmed đã có)"
        return
    fi

    echo -n "Đánh giá chất lượng trình tự (Q1, Q2, Q3): "
    read a

    case "$a" in
        Q1)
            log_run "fastp – Q1 (chất lượng tốt, không trim mạnh)"
            conda run --no-capture-output -n pipeline_env fastp \
              -i SRread_1.fastq.gz -o SRread_1_trimmed.fastq.gz \
              -I SRread_2.fastq.gz -O SRread_2_trimmed.fastq.gz \
              -j QC_report/report.json \
              -h QC_report/report.html \
              -w 8
            ;;
        Q2)
            log_run "fastp – Q2 (chất lượng trung bình)"
            conda run --no-capture-output -n pipeline_env fastp \
              -i SRread_1.fastq.gz -o SRread_1_trimmed.fastq.gz \
              -I SRread_2.fastq.gz -O SRread_2_trimmed.fastq.gz \
              --cut_tail --cut_mean_quality 25 --cut_window_size 4 \
              --length_required 100 \
              -j QC_report/report.json -h QC_report/report.html \
              -w 8
            ;;
        Q3)
            log_run "fastp – Q3 (chất lượng thấp)"
            conda run --no-capture-output -n pipeline_env fastp \
              -i SRread_1.fastq.gz -o SRread_1_trimmed.fastq.gz \
              -I SRread_2.fastq.gz -O SRread_2_trimmed.fastq.gz \
              --cut_right --cut_mean_quality 20 --cut_window_size 4 \
              --average_qual 25 --qualified_quality_phred 15 --unqualified_percent_limit 30 \
              --length_required 100 \
              -j QC_report/report.json -h QC_report/report.html \
              -w 8
            ;;
        *)
            log_err "Lựa chọn không hợp lệ: '$a'. Chỉ chấp nhận Q1, Q2, Q3."
            exit 1
            ;;
    esac

    log_done "Trimming xong. File: SRread_1_trimmed.fastq.gz, SRread_2_trimmed.fastq.gz, QC_report/report.html"
}

# ══════════════════════════════════════════════════════════
# 3: FASTQC CHO READ ĐÃ TRIM
# ══════════════════════════════════════════════════════════
step_qc2() {
    if exists "QC_report/SRread_1_trimmed_fastqc.html" && exists "QC_report/SRread_2_trimmed_fastqc.html"; then
        log_skip "FastQC read đã trim (report đã có)"
    else
        log_run "FastQC cho read đã trim"
        conda run --no-capture-output -n pipeline_env fastqc -o QC_report/ SRread_1_trimmed.fastq.gz
        conda run --no-capture-output -n pipeline_env fastqc -o QC_report/ SRread_2_trimmed.fastq.gz
    fi
    log_done "FastQC read đã trim xong. File: QC_report/SRread_1_trimmed_fastqc.html, QC_report/SRread_2_trimmed_fastqc.html"
}

# ══════════════════════════════════════════════════════════
# 4: ASSEMBLY VỚI UNICYCLER (chuyển file rồi chạy, log hiển thị trực tiếp)
# ══════════════════════════════════════════════════════════
step_assembly() {
    if exists "assembly_fastq/SRread_1_for_assembly.fastq.gz" && exists "assembly_fastq/SRread_2_for_assembly.fastq.gz"; then
        log_skip "Chuyển file vào assembly_fastq/ (đã có)"
    else
        log_run "Chuyển trimmed read vào assembly_fastq/"
        mv SRread_1_trimmed.fastq.gz assembly_fastq/SRread_1_for_assembly.fastq.gz
        mv SRread_2_trimmed.fastq.gz assembly_fastq/SRread_2_for_assembly.fastq.gz
    fi

    if exists "assembly/assembly.fasta"; then
        log_skip "Unicycler assembly (assembly/assembly.fasta đã có)"
    else
        log_run "Assembly với Unicycler (log chạy hiển thị trực tiếp bên dưới)"
        conda run --no-capture-output -n unicycler_env \
          unicycler \
            -1 assembly_fastq/SRread_1_for_assembly.fastq.gz \
            -2 assembly_fastq/SRread_2_for_assembly.fastq.gz \
            --min_fasta_length 150 \
            --verbosity 2 \
            -o assembly
    fi

    log_done "Assembly xong. File: assembly/assembly.fasta (thư mục assembly/)"
}

# ══════════════════════════════════════════════════════════
# 5: BUSCO
# ══════════════════════════════════════════════════════════
step_busco() {
    if ls busco_result/short_summary*.txt 2>/dev/null | grep -q .; then
        log_skip "BUSCO (short_summary đã có trong busco_result/)"
    else
        log_run "BUSCO đánh giá độ đầy đủ genome (log chạy hiển thị trực tiếp bên dưới)"
        conda run --no-capture-output -n busco_env \
          busco \
            -i assembly/assembly.fasta \
            -o busco_result \
            -l bacteria_odb10 \
            -m genome \
            --download_path busco_download \
            -c 8 \
            --force
    fi
    log_done "BUSCO xong. File: busco_result/short_summary*.txt (thư mục busco_result/)"
}

# ══════════════════════════════════════════════════════════
# 6: QUAST
# ══════════════════════════════════════════════════════════
step_quast() {
    if exists "quast_result/report.html"; then
        log_skip "QUAST (report.html đã có trong quast_result/)"
    else
        log_run "QUAST thống kê chất lượng assembly (log chạy hiển thị trực tiếp bên dưới)"
        conda run --no-capture-output -n pipeline_env quast \
          -o quast_result \
          -t 8 \
          assembly/assembly.fasta
    fi
    log_done "QUAST xong. File: quast_result/report.html (thư mục quast_result/)"
}

run_all() {
    step_install
    step_rename
    step_qc1
    step_trim
    step_qc2
    step_assembly
    step_busco
    step_quast
}

# ══════════════════════════════════════════════════════════
# MENU
# ══════════════════════════════════════════════════════════
show_menu() {
    echo ""
    echo "========================================"
    echo "   PIPELINE ASSEMBLY - CHỌN CHỨC NĂNG"
    echo "========================================"
    echo "  -1  Cài đặt môi trường & công cụ"
    echo "   0  Đổi tên file fastq.gz"
    echo "   1  FastQC - đọc gốc"
    echo "   2  Trimming (fastp)"
    echo "   3  FastQC - đọc sau trim"
    echo "   4  Assembly (Unicycler)"
    echo "   5  BUSCO"
    echo "   6  QUAST"
    echo "   a  Chạy tất cả các bước"
    echo "   q  Thoát"
    echo "========================================"
}

run_choice() {
    case "$1" in
        -1) step_install ;;
        0)  step_rename ;;
        1)  step_qc1 ;;
        2)  step_trim ;;
        3)  step_qc2 ;;
        4)  step_assembly ;;
        5)  step_busco ;;
        6)  step_quast ;;
        a|A) run_all ;;
        q|Q) exit 0 ;;
        *)
            log_err "Lựa chọn không hợp lệ: '$1'"
            exit 1
            ;;
    esac
}

# ══════════════════════════════════════════════════════════
# MAIN
# ══════════════════════════════════════════════════════════
run_and_report() {
    local choice="$1"
    local t0 t1 dur
    t0=$(date +%s)

    init_dirs
    run_choice "$choice"

    t1=$(date +%s)
    dur=$((t1 - t0))
    echo ""
    echo "----------------------------------------"
    echo "KẾT THÚC: $(date '+%Y-%m-%d %H:%M:%S')"
    echo "Thời gian chạy: ${dur} giây."
    echo "----------------------------------------"
}

CHOICE="$1"

if [ -n "$CHOICE" ]; then
    # Có truyền tham số dòng lệnh -> chạy 1 lần rồi thoát (không hiện menu)
    run_and_report "$CHOICE"
else
    # Không truyền tham số -> hiện menu, chạy xong quay lại menu, đến khi chọn q
    while true; do
        show_menu
        echo -n "Nhập lựa chọn: "
        read CHOICE
        run_and_report "$CHOICE"
        echo ""
    done
fi
