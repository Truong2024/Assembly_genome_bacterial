#!/usr/bin/env bash
set -uo pipefail

# ==============================================================
# PIPELINE: Bakta - GTDB-TK - Pan-genome - Sàng lọc AMR/Yếu tố độc lực
# ==============================================================
#   Phần -1: Cài đặt & kiểm tra công cụ (conda env, binary, CSDL)
#   Phần 0: Kiểm tra file assembly.fasta
#   Phần 1: Chú giải gen với Bakta
#   Phần 2: Phân loại GTDB-TK
#   Phần 3: Tải genome tham chiếu (accession) cho pan-genome
#   Phần 4: Pan-genome (Prokka + Roary + IQ-TREE + biểu đồ)
#   Phần 5: Sàng lọc AMR / Yếu tố độc lực (Abricate)
# ==============================================================

# ---------------- Cấu hình chung (chỉnh trước khi chạy) ----------------
GENOME="${GENOME:-genome/assembly.fasta}"          # genome của mình
STRAIN="${STRAIN:-strain1}"                        # tên chủng, dùng làm tên file
GENUS="${GENUS:-Bacillus}"                         # chi (dùng cho Prokka/Bakta)
SPECIES="${SPECIES:-subtilis}"                     # loài (dùng cho Prokka/Bakta)
THREADS="${THREADS:-6}"                            # số luồng CPU
GTDBTK_DATA_PATH_DEFAULT="/mnt/e/UniGenLab/DATABASE_ALL/release232"

# ---------------- Cấu hình riêng cho Bakta (chỉnh để tối ưu chú giải) ----------------
BAKTA_DB_PATH_DEFAULT="/mnt/e/UniGenLab/DATABASE_ALL/bakta_DB/db"
BAKTA_COMPLETE="${BAKTA_COMPLETE:-false}"          # genome đang là draft (nhiều contigs) -> không bật --complete
BAKTA_COMPLIANT="${BAKTA_COMPLIANT:-true}"         # cần nộp NCBI/ENA -> luôn bật --compliant + locus-tag
BAKTA_GRAM="${BAKTA_GRAM:-?}"                      # để mặc định "?" -> script sẽ hỏi + hoặc - khi chạy Phần 1
BAKTA_LOCUS_TAG="${BAKTA_LOCUS_TAG:-$STRAIN}"      # locus tag prefix (tự động viết hoa khi dùng --compliant)

# ---------------- Hàm tiện ích ----------------
log_run()  { echo -e "\033[1;34m[Running...]\033[0m   $*"; }
log_skip() { echo -e "\033[1;33m[Skipped]\033[0m $*"; }
log_err()  { echo -e "\033[1;31m[Error]\033[0m   $*" >&2; }
log_done() { echo -e "\033[1;32m[Done]\033[0m    $*"; }
exists()   { [ -e "$1" ]; }

mkdir -p genome gtdbtk bakta_output pangenome/genomes pangenome/gff_files pangenome/prokka screening_output

# ==============================================================
# PHẦN -1 - Cài đặt & kiểm tra công cụ
# ==============================================================
# Kiểm tra từng môi trường conda + lệnh cần dùng trong đó. Nếu thiếu,
# hỏi người dùng có muốn cài luôn không (tên env/tool giữ nguyên như
# phần còn lại của script, không đổi tên).
check_or_install_tool() {
    local env_name="$1" tool_cmd="$2" pkg="$3"
    echo -n "  - Môi trường '$env_name' (cần lệnh: $tool_cmd) ... "

    if ! conda env list | awk '{print $1}' | grep -qx "$env_name"; then
        echo -e "\033[1;31mCHƯA CÓ\033[0m"
        read -rp "    Tạo môi trường '$env_name' và cài '$pkg' ngay? (y/N): " ans
        if [[ "$ans" =~ ^[Yy]$ ]]; then
            conda create -y -n "$env_name" -c conda-forge -c bioconda "$pkg"
        else
            log_err "    Bỏ qua. Các phần dùng '$env_name' sẽ lỗi cho tới khi được cài."
            return 1
        fi
        return 0
    fi

    if conda run -n "$env_name" bash -lc "command -v $tool_cmd" >/dev/null 2>&1; then
        echo -e "\033[1;32mOK\033[0m"
        return 0
    fi

    echo -e "\033[1;33mmôi trường có nhưng thiếu lệnh '$tool_cmd'\033[0m"
    read -rp "    Cài '$pkg' vào môi trường '$env_name' ngay? (y/N): " ans
    if [[ "$ans" =~ ^[Yy]$ ]]; then
        conda install -y -n "$env_name" -c conda-forge -c bioconda "$pkg"
    else
        log_err "    Bỏ qua. '$tool_cmd' sẽ không chạy được trong '$env_name'."
        return 1
    fi
}

partm1_check_tools() {
    if ! command -v conda >/dev/null 2>&1; then
        log_err "Phần -1: không tìm thấy lệnh 'conda' trong PATH. Cần cài Miniconda/Anaconda trước."
        return 1
    fi

    log_run "Phần -1: kiểm tra các môi trường conda & công cụ dòng lệnh"
    check_or_install_tool "bakta_env" "bakta"    "bioconda::bakta"
    check_or_install_tool "gtdbtk"    "gtdbtk"   "bioconda::gtdbtk"
    check_or_install_tool "taxonomy"  "datasets" "conda-forge::ncbi-datasets-cli"
    check_or_install_tool "annotate"  "prokka"   "bioconda::prokka"
    check_or_install_tool "pangenome" "roary"    "bioconda::roary"
    check_or_install_tool "pangenome" "iqtree2"  "bioconda::iqtree"
    check_or_install_tool "screening" "abricate" "bioconda::abricate"

    echo -n "  - Công cụ hệ thống 'unzip' ... "
    command -v unzip >/dev/null 2>&1 && echo -e "\033[1;32mOK\033[0m" || echo -e "\033[1;31mTHIẾU\033[0m (cần cho Phần 3)"

    echo -n "  - Công cụ hệ thống 'python3' ... "
    command -v python3 >/dev/null 2>&1 && echo -e "\033[1;32mOK\033[0m" || echo -e "\033[1;31mTHIẾU\033[0m (cần cho Phần 4.5)"

    echo -n "  - Script 'roary_plots.py' (thư mục hiện tại) ... "
    if exists "roary_plots.py"; then
        echo -e "\033[1;32mOK\033[0m"
    else
        echo -e "\033[1;33mKHÔNG THẤY\033[0m (copy từ mã nguồn Roary, thư mục contrib/roary_plots/, vào thư mục chạy script nếu muốn có Phần 4.5)"
    fi

    echo -n "  - CSDL Bakta tại $BAKTA_DB_PATH_DEFAULT ... "
    if [ -d "$BAKTA_DB_PATH_DEFAULT" ] && [ -n "$(ls -A "$BAKTA_DB_PATH_DEFAULT" 2>/dev/null)" ]; then
        echo -e "\033[1;32mOK\033[0m"
    else
        echo -e "\033[1;33mCHƯA THẤY\033[0m (tải bằng lệnh 'bakta_db download' trong môi trường bakta_env trước khi chạy Phần 1)"
    fi

    echo -n "  - CSDL GTDB-TK tại $GTDBTK_DATA_PATH_DEFAULT ... "
    if [ -d "$GTDBTK_DATA_PATH_DEFAULT" ] && [ -n "$(ls -A "$GTDBTK_DATA_PATH_DEFAULT" 2>/dev/null)" ]; then
        echo -e "\033[1;32mOK\033[0m"
    else
        echo -e "\033[1;33mCHƯA THẤY\033[0m (Phần 2 sẽ tự động tải bằng download-db.sh khi chạy trong môi trường gtdbtk)"
    fi

    log_done "Phần -1: kiểm tra công cụ hoàn tất (xem chi tiết OK/THIẾU ở trên)."
}

# ==============================================================
# PHẦN 0 - Kiểm tra file assembly.fasta
# ==============================================================
part0_check_assembly() {
    if exists "$GENOME"; then
        log_skip "Phần 0: đã tìm thấy assembly ($GENOME)"
        return 0
    else
        log_err "Phần 0: không tìm thấy $GENOME. Kiểm tra lại biến GENOME."
        return 1
    fi
}

# ==============================================================
# PHẦN 1 - Chú giải gen với Bakta
# ==============================================================
part1_bakta() {
    part0_check_assembly || return 1

    if exists "bakta_output/${STRAIN}.tsv"; then
        log_skip "Phần 1: Bakta (đã chạy trước đó)"
        log_done "Phần 1: kết quả đã có sẵn tại bakta_output/ (prefix: ${STRAIN})"
        return 0
    fi

    # Hỏi Gram (+/-) trực tiếp nếu chưa được set sẵn qua biến môi trường BAKTA_GRAM
    while [ "$BAKTA_GRAM" != "+" ] && [ "$BAKTA_GRAM" != "-" ]; do
        echo -n "Nhập Gram của chủng (chỉ nhập + hoặc -): "
        read -r BAKTA_GRAM
    done

    export BAKTA_DB="${BAKTA_DB:-$BAKTA_DB_PATH_DEFAULT}"

    local bakta_args=(
        --db "$BAKTA_DB"
        --output bakta_output
        --prefix "$STRAIN"
        --genus "$GENUS"
        --species "$SPECIES"
        --strain "$STRAIN"
        --threads "$THREADS"
        --gram "$BAKTA_GRAM"
        --force
    )

    [ "$BAKTA_COMPLETE" = "true" ] && bakta_args+=(--complete)
    if [ "$BAKTA_COMPLIANT" = "true" ]; then
        bakta_args+=(--compliant --locus-tag "${BAKTA_LOCUS_TAG^^}")
    fi

    log_run "Phần 1: chú giải gen bằng Bakta (DB: $BAKTA_DB, Gram: $BAKTA_GRAM, compliant: $BAKTA_COMPLIANT)"
    conda run -n bakta_env bakta "${bakta_args[@]}" "$GENOME"
    local rc=$?
    if [ $rc -eq 0 ]; then
        log_done "Phần 1: Chú giải gen hoàn tất. Kết quả lưu tại: bakta_output/ (prefix: ${STRAIN})"
    else
        log_err "Phần 1: Bakta kết thúc với lỗi (exit code $rc)."
    fi
    return $rc
}

# ==============================================================
# PHẦN 2 - Phân loại GTDB-TK
# ==============================================================
part2_gtdbtk() {
    part0_check_assembly || return 1

    if exists "gtdbtk/classify/gtdbtk.bac120.summary.tsv" || exists "gtdbtk/gtdbtk.bac120.summary.tsv"; then
        log_skip "Phần 2: GTDB-TK classify_wf (đã chạy trước đó)"
        log_done "Phần 2: kết quả đã có sẵn tại gtdbtk/"
        return 0
    fi

    conda activate gtdbtk

    export GTDBTK_DATA_PATH="${GTDBTK_DATA_PATH:-$GTDBTK_DATA_PATH_DEFAULT}"
    echo "GTDBTK_DATA_PATH=$GTDBTK_DATA_PATH"

    if [ ! -d "$GTDBTK_DATA_PATH" ] || [ -z "$(ls -A "$GTDBTK_DATA_PATH" 2>/dev/null)" ]; then
        log_run "Phần 2.1: tải CSDL GTDB (release232)"
        download-db.sh
    else
        log_skip "Phần 2.1: CSDL GTDB đã tồn tại tại $GTDBTK_DATA_PATH"
    fi

    log_run "Phần 2.2: chạy gtdbtk classify_wf"
    conda run -n gtdbtk gtdbtk classify_wf --genome_dir genome --out_dir gtdbtk -x fasta --cpus "$THREADS" \
        --scratch_dir gtdbtk/gtdbtk_scratch
    local rc=$?
    if [ $rc -eq 0 ]; then
        log_done "Phần 2: Phân loại GTDB-TK hoàn tất. Kết quả lưu tại: gtdbtk/ (xem gtdbtk/classify/gtdbtk.bac120.summary.tsv hoặc gtdbtk/gtdbtk.bac120.summary.tsv)"
    else
        log_err "Phần 2: GTDB-TK kết thúc với lỗi (exit code $rc)."
    fi
    return $rc
}

# ==============================================================
# PHẦN 3 - Tải genome tham chiếu (accession) cho pan-genome
# ==============================================================
# Tách riêng khỏi Phần 4 (phân tích Roary) để có thể tải trước, kiểm
# tra dữ liệu, rồi mới chạy phân tích mà không phải tải lại.
part3_download_refs() {
    part0_check_assembly || return 1

    if exists "pangenome/genomes_downloaded.flag"; then
        log_skip "Phần 3: tải genome tham chiếu (đã tải trước đó)"
        log_done "Phần 3: dữ liệu đã có sẵn tại bsubtilis_refs/"
        return 0
    fi

    echo -n "Nhập danh sách accession genome tham chiếu (cách nhau bởi dấu phẩy, vd: GCF_000009045.1,GCF_000009048.1): "
    read -r acc_list

    if [ -z "$acc_list" ]; then
        log_err "Chưa nhập accession nào. Hủy Phần 3."
        return 1
    fi

    log_run "Phần 3: tải genome accession: $acc_list"
    conda run -n taxonomy datasets download genome accession "$acc_list" \
        --include genome --filename bsubtilis_refs.zip
    local rc=$?
    if [ $rc -ne 0 ]; then
        log_err "Phần 3: tải genome thất bại (exit code $rc)."
        return $rc
    fi

    unzip -o bsubtilis_refs.zip -d bsubtilis_refs/
    touch pangenome/genomes_downloaded.flag
    log_done "Phần 3: Tải & giải nén genome tham chiếu hoàn tất. Kết quả lưu tại: bsubtilis_refs/ (đánh dấu hoàn tất: pangenome/genomes_downloaded.flag)"
}

# ==============================================================
# PHẦN 4 - Pan-genome (Prokka + Roary + IQ-TREE + biểu đồ)
# ==============================================================
part4_pangenome_analysis() {
    part0_check_assembly || return 1

    if ! exists "pangenome/genomes_downloaded.flag"; then
        log_err "Phần 4: chưa có genome tham chiếu. Hãy chạy Phần 3 trước."
        return 1
    fi

    # 4.1 Chuẩn bị thư mục genome (genome của mình + genome tham chiếu)
    if ls pangenome/genomes/*.fasta >/dev/null 2>&1 || ls pangenome/genomes/*.fna >/dev/null 2>&1; then
        log_skip "Phần 4.1: copy genome vào pangenome/genomes/"
    else
        log_run "Phần 4.1: copy genome vào pangenome/genomes/"
        cp "$GENOME" pangenome/genomes/"${STRAIN}".fasta
        cp bsubtilis_refs/ncbi_dataset/data/*/*.fna pangenome/genomes/ 2>/dev/null || true
    fi

    # 4.2 Prokka đồng nhất cho tất cả genome
    if ls pangenome/gff_files/*.gff >/dev/null 2>&1; then
        log_skip "Phần 4.2: Prokka đồng nhất cho pan-genome"
    else
        log_run "Phần 4.2: Prokka đồng nhất cho từng genome trong pangenome/genomes/"
        for f in pangenome/genomes/*.fasta pangenome/genomes/*.fna; do
            [ -e "$f" ] || continue
            name=$(basename "$f" | sed 's/\.[^.]*$//')
            if exists "pangenome/prokka/$name/$name.gff"; then
                log_skip "  - $name (đã chú giải)"
            else
                log_run "  - Prokka cho $name"
                conda run -n annotate prokka --outdir pangenome/prokka/"$name" --prefix "$name" \
                    --kingdom Bacteria --genus "$GENUS" --species "$SPECIES" \
                    --cpus "$THREADS" --force "$f"
            fi
            cp pangenome/prokka/"$name"/"$name".gff pangenome/gff_files/ 2>/dev/null || true
        done
    fi

    # 4.3 Roary
    if exists "pangenome/roary_output/summary_statistics.txt"; then
        log_skip "Phần 4.3: Roary"
    else
        log_run "Phần 4.3: Roary pan-genome analysis"
        rm -rf pangenome/roary_output  # roary yêu cầu thư mục output chưa tồn tại
        conda run -n pangenome roary -e --mafft -p "$THREADS" \
            -f pangenome/roary_output pangenome/gff_files/*.gff
        if [ $? -ne 0 ]; then
            log_err "Phần 4.3: Roary kết thúc với lỗi."
            return 1
        fi
    fi

    # 4.4 Cây phát sinh loài (IQ-TREE)
    if exists "pangenome/core_tree.treefile"; then
        log_skip "Phần 4.4: cây phát sinh loài (IQ-TREE)"
    else
        log_run "Phần 4.4: IQ-TREE từ core gene alignment"
        conda run -n pangenome iqtree2 \
            -s pangenome/roary_output/core_gene_alignment.aln \
            -m GTR+G -bb 1000 -nt "$THREADS" -pre pangenome/core_tree
        if [ $? -ne 0 ]; then
            log_err "Phần 4.4: IQ-TREE kết thúc với lỗi."
            return 1
        fi
    fi

    # 4.5 Biểu đồ pan-genome (roary_plots.py)
    if exists "pangenome/roary_plots"; then
        log_skip "Phần 4.5: biểu đồ Roary"
    else
        if exists "roary_plots.py"; then
            log_run "Phần 4.5: vẽ biểu đồ Roary (pie chart, presence/absence)"
            python3 roary_plots.py \
                pangenome/roary_output/accessory_binary_genes.fa.newick \
                pangenome/roary_output/gene_presence_absence.csv \
                --outdir pangenome/roary_plots
        else
            log_err "Thiếu roary_plots.py trong thư mục hiện tại — bỏ qua Phần 4.5."
        fi
    fi

    log_done "Phần 4: Pan-genome hoàn tất. Kết quả lưu tại: pangenome/roary_output/ (thống kê + gene_presence_absence.csv), pangenome/core_tree.* (cây IQ-TREE), pangenome/roary_plots/ (biểu đồ, nếu có roary_plots.py)"
}

# ==============================================================
# PHẦN 5 - Sàng lọc AMR / Yếu tố độc lực (Abricate)
# ==============================================================
part5_screening() {
    part0_check_assembly || return 1

    if exists "screening_output/amr_resfinder.tsv"; then
        log_skip "Phần 5.1: AMR (ResFinder)"
    else
        log_run "Phần 5.1: AMR (ResFinder)"
        conda run -n screening abricate --db resfinder "$GENOME" > screening_output/amr_resfinder.tsv
    fi

    if exists "screening_output/amr_card.tsv"; then
        log_skip "Phần 5.2: AMR (CARD)"
    else
        log_run "Phần 5.2: AMR (CARD)"
        conda run -n screening abricate --db card "$GENOME" > screening_output/amr_card.tsv
    fi

    if exists "screening_output/virulence_vfdb.tsv"; then
        log_skip "Phần 5.3: Yếu tố độc lực (VFDB)"
    else
        log_run "Phần 5.3: Yếu tố độc lực (VFDB)"
        conda run -n screening abricate --db vfdb "$GENOME" > screening_output/virulence_vfdb.tsv
    fi

    if exists "screening_output/amr_summary.tsv"; then
        log_skip "Phần 5.4: tổng hợp AMR summary"
    else
        log_run "Phần 5.4: tổng hợp AMR summary"
        conda run -n screening abricate --summary screening_output/amr_resfinder.tsv > screening_output/amr_summary.tsv
    fi

    log_done "Phần 5: Sàng lọc AMR/Yếu tố độc lực hoàn tất. Kết quả lưu tại: screening_output/ (amr_resfinder.tsv, amr_card.tsv, virulence_vfdb.tsv, amr_summary.tsv)"
}

# ==============================================================
# MENU
# ==============================================================
run_choice() {
    case "$1" in
        -1) partm1_check_tools ;;
        0) part0_check_assembly ;;
        1) part1_bakta ;;
        2) part2_gtdbtk ;;
        3) part3_download_refs ;;
        4) part4_pangenome_analysis ;;
        5) part5_screening ;;
        6)
            part0_check_assembly && part1_bakta && part2_gtdbtk && part3_download_refs && part4_pangenome_analysis && part5_screening
            ;;
        *) log_err "Lựa chọn không hợp lệ: $1" ;;
    esac
}

show_menu() {
    echo ""
    echo "================================================================"
    echo "   MENU - BAKTA / GTDB-TK / PAN-GENOME / SÀNG LỌC AMR"
    echo "================================================================"
    echo " -1) Cài đặt & kiểm tra công cụ (conda env, binary, CSDL)"
    echo "  0) Kiểm tra file assembly.fasta"
    echo "  1) Chú giải gen với Bakta"
    echo "  2) Phân loại GTDB-TK"
    echo "  3) Tải genome tham chiếu (accession) cho pan-genome"
    echo "  4) Pan-genome (Prokka + Roary + IQ-TREE + biểu đồ)"
    echo "  5) Sàng lọc AMR / Yếu tố độc lực (Abricate)"
    echo "  6) Chạy tất cả (0 -> 5)"
    echo "  7) Thoát"
    echo "----------------------------------------------------------------"
    echo "  Có thể chọn nhiều mục cùng lúc, cách nhau bởi khoảng trắng."
    echo "  Vd: 1 2 3"
    echo "  Gợi ý: chạy mục -1 trước khi dùng pipeline lần đầu."
}

main() {
    while true; do
        show_menu
        read -rp "Chọn chức năng: " -a choices

        [ "${#choices[@]}" -eq 0 ] && continue

        for c in "${choices[@]}"; do
            if [ "$c" = "7" ]; then
                echo "Thoát."
                exit 0
            fi
            run_choice "$c"
        done
    done
}

main "$@"
