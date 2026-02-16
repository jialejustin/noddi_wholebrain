#!/usr/bin/env python3
"""
Produces whole-brain NODDI metrics from noddi_reg ROI outputs

Usage:
    noddi_wholebrain [options] --participants_tsv <participants_tsv> --noddi_reg_dir <noddi_reg_dir> --output_dir <output_dir>

Arguments:
    --participants_tsv <participants_tsv>    Participants.tsv file
    --noddi_reg_dir <noddi_reg_dir>          Input directory from noddi_reg
    --output_dir <output_dir>                Output directory

Options:
    --debug                          Debug logging
    -h, --help                       Prints this message

DETAILS:

Extracts whole-brain NODDI metrics from noddi_reg ROI outputs for a list of participants.
The whole-brain metric is calculated as an average of the mean NODDI metric in each ROI weighted by it's masked voxels.

"""
from docopt import docopt

import os.path

from bids import BIDSLayout
import pandas as pd

import logging

logger = logging.getLogger(os.path.basename(__file__))

# read freesurfer tissue types template globally
freesurfer_tissue_types = pd.read_csv(
    "/scratch/juyu/ScanD_SCC/templates/desc-FreeSurferAll_dseg_with_tissue_type.tsv",
    sep="\t",
)[["name", "tissue_type"]]


# calculate whole-brain NODDI
def get_wholebrain_noddi(parc_files):
    wb_noddi_list = []
    for parc_file in parc_files:
        parc_df = parc_file.get_df()
        desc = parc_file.get_entities()["desc"]

        # filter to GM if freesurfer
        if desc == "aparcaseg":
            parc_df = parc_df.merge(freesurfer_tissue_types, how="left")
            parc_df = parc_df[parc_df["tissue_type"] == "GM"]

        # get whole-brain GM NODDI via weighted average of all GM ROIs by voxel count

        total_vx = parc_df["n_vx_masked"].sum()

        metrics = {
            "parcellation": desc,
            "whole_od": (parc_df["od_mean"] * parc_df["n_vx_masked"]).sum() / total_vx,
            "whole_icvf": (parc_df["icvf_mean"] * parc_df["n_vx_masked"]).sum() / total_vx,
            "whole_isovf": (parc_df["isovf_mean"] * parc_df["n_vx_masked"]).sum() / total_vx,
        }
        wb_noddi_list.append(metrics)

    return wb_noddi_list


def main():
    # handle docopt arguments
    arguments = docopt(__doc__)

    logger.setLevel(logging.WARNING)

    if arguments["--debug"]:
        logger.setLevel(logging.DEBUG)

    logger.info(arguments)

    participants_tsv = arguments["--participants_tsv"]
    noddi_reg_dir = arguments["--noddi_reg_dir"]
    output_dir = arguments["--output_dir"]


    participants_list = pd.read_csv(participants_tsv, sep="\t")["participant_id"]
    noddireg_layout = BIDSLayout(noddi_reg_dir, is_derivative=True, validate=False)

    all_participant_results = []
    for participant in participants_list:
        # strip 'sub-' from tsv codes
        subject = participant.replace("sub-", "")
        parc_files = noddireg_layout.get(
            subject=subject,
            datatype="dwi",
            model="noddi",
            suffix="results",
            extension=".tsv",
        )

        if not parc_files:
            logger.error(f"No noddi_reg results found for {participant}")
            continue

        # retrieve list of dicts with whole-brain noddi for single subject
        wb_noddi_list = get_wholebrain_noddi(parc_files)

        # add participant to results for later
        for wb_noddi_metric in wb_noddi_list:
            wb_noddi_metric["participant_id"] = participant
            all_participant_results.append(wb_noddi_metric)

    # create big dataframe with all results
    df_all = pd.DataFrame(all_participant_results)

    # split by parcellation and save to csv
    for parc_type, parc_df in df_all.groupby("parcellation"):
        output_path = os.path.join(output_dir, f"desc-{parc_type}_wholebrainnoddi.csv")
        # remove redundant parcellation column
        parc_df = parc_df.drop(columns=["parcellation"])
        # reorder and save
        parc_df[["participant_id", "whole_od", "whole_icvf", "whole_isovf"]].to_csv(
            output_path, index=False
        )


if __name__ == "__main__":
    main()
