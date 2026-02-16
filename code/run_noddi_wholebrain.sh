BASEDIR=${PWD}

PARTICIPANTS_TSV=${BASEDIR}/data/local/bids/participants.tsv
NODDI_REG_DIR=${BASEDIR}/data/local/derivatives/noddi_reg
OUTPUT_DIR=${BASEDIR}/data/local/derivatives/noddi_wholebrain

mkdir -p ${OUTPUT_DIR}

python code/noddi_wholebrain.py \
    --participants_tsv "${PARTICIPANTS_TSV}" \
    --noddi_reg_dir "${NODDI_REG_DIR}" \
    --output_dir "${OUTPUT_DIR}"