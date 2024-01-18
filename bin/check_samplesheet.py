#!/usr/bin/env python

import sys
import argparse
import re
import pandas as pd

def parse_args(args=None):
    Description = "Reformat samplesheet file and check its contents."

    Epilog = "Example usage: python check_samplesheet.py <FILE_IN> <FILE_OUT>"

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument("FILE_IN", help="Input samplesheet file.")
    parser.add_argument("FILE_OUT", help="Output file.")
    return parser.parse_args(args)

def print_error(error, context="Line", context_str=""):
    error_str = "ERROR: Please check samplesheet -> {}".format(error)
    if context != "" and context_str != "":
        error_str = "ERROR: Please check samplesheet -> {}\n{}: '{}'".format(
            error, context.strip(), context_str.strip()
        )
    print(error_str)
    sys.exit(1)

def match_legal_pattern(label):
    '''Helper function to check if labels are legal'''
    legal_pattern = r"^[a-zA-Z][a-zA-Z0-9_]*$"
    if not re.match(legal_pattern, label):
        print_error("Sample/Group label contain illegal characters or does not start with letters", "Label", label)
    # Mostly MultiQC reserved strings, with some reserved by the pipeline
    illegal_patterns = ["_tophat", "ReadsPerGene", "_star_aligned", "_fastqc", "_counts", "Aligned", "_slamdunk", 
                        "_bismark", "_SummaryStatistics", "_duprate", "_vep", "ccs", "_NanoStats", 
                        r"_trimmed$", r"_val$", r"_mqc$", r"short_summary_$", r"^short_summary_", r"_summary$", r"_matrix$",
                        r"_$", r"_R1$", r"_R2$"]
    for p in illegal_patterns:
        if re.match(p, label):
            print_error("Sample/Group label {} contain string reserved by pipeline or MultiQC. This may cause a problem.".format(label),
                        "Reserved string", p)
    return True

def check_fastq_suffix(filename):
    '''Helper function to check if filenames have the correct suffix'''
    FQ_EXTENSIONS = (".fq.gz", ".fastq.gz")
    if not filename.endswith(FQ_EXTENSIONS):
        print_error("FASTQ path has invalid extension", "Path", filename)
    return True

def check_all_se_or_all_pe(group):
    '''Helper function to check if all runs of same sample are either all single-ended or all paired-ended'''
    return group['read_2'].count() == 0 or group['read_2'].count() == group['read_1'].count()

def check_samplesheet(file_in, file_out):
    """
    This function checks that the samplesheet follows the following structure:

    sample,read_1,read_2,group,run_accession
    sample1,s1_run1_R1.fastq.gz,s1_run1_R2.fastq.gz,groupA,run1
    sample1,s1_run2_R1.fastq.gz,s1_run2_R2.fastq.gz,groupA,run2
    sample2,s2_run1_R1.fastq.gz,,groupB,,
    sample3,s3_run1_R1.fastq.gz,s3_run1_R2.fastq.gz,groupB,,

    run_accession only required for rows with duplicated sample names indicating those reads should be combined
    """

    design = pd.read_csv(file_in, index_col=False)

    # Check column names
    HEADER = {"sample", "read_1", "read_2", "group"}
    assert HEADER.issubset(set(design.columns)), "Design file must contain {} columns".format(','.join(HEADER))

    # Check if all sample labels legal
    assert design["sample"].map(match_legal_pattern).all()

    # Make sure group labels are either all missing or all present
    assert design["group"].isna().all() or design["group"].notnull().all(), "Group labels missing in some samples but not others!"

    # Check if all group labels legal
    assert design["group"].map(match_legal_pattern, na_action="ignore").all()

    # Check FASTQ locations
    assert design["read_1"].all(), "Read 1 path cannot be missing!"

    # Check FASTQ extension
    assert design["read_1"].map(check_fastq_suffix).all()
    assert design["read_2"].map(check_fastq_suffix, na_action="ignore").all()
    
    # Make sure there is no duplication of FASTQ locations
    assert design["read_1"].duplicated().any() == False, "There are duplications within Read 1 paths!"
    assert design["read_2"].dropna().duplicated().any() == False, "There are duplications within Read 2 paths!"
    
    # Make sure group labels are consistent for the same sample, different runs
    group_label_consistent = design.groupby("sample")["group"].apply(lambda x:len(x.unique())==1)
    if not group_label_consistent.all():
        inconsistent_group_labels = group_label_consistent[~group_label_consistent].index.tolist()
        print_error("Group labels for samples not consistent across different runs", "Samples", ",".join(inconsistent_group_labels))

    # Make sure there are no mixtures of single and paired-end data within the same sample, different runs
    sample_all_se_or_pe = design.groupby("sample").apply(check_all_se_or_all_pe)
    if not sample_all_se_or_pe.all():
        bad_samples = sample_all_se_or_pe[~sample_all_se_or_pe].index.tolist()
        print_error("Single-end and paired-end data cannot be mixed for the same sample", "Samples", ",".join(bad_samples))

    # Make sure run_accession is present when there are duplicated samples
    if design["sample"].duplicated().any():
        if "run_accession" not in design.columns:
            duplicated_samples = design["sample"][design["sample"].duplicated()].tolist()
            print_error("run_accesssion column must exist when there are duplicated sample labels", "Samples", ",".join(duplicated_samples))

        # Make sure run_accession legal
        assert design["run_accession"].map(match_legal_pattern, na_action="ignore").all()
    
        # Make sure run_accession is not duplicated for the same sample
        accession_counts = design.groupby("sample")["run_accession"].agg(["size","nunique"])
        for idx, row in accession_counts.iterrows():
            if row["size"]>1 and row["size"] != row["nunique"]:
                print_error("run_accession missing or not unique for the same sample", "Sample", idx)

    ###################################################
    # PREPARE for OUTPUT
    # Add a "run_accession" column if one doesn't exist
    if "run_accession" not in design.columns:
        design["run_accession"] = ""
    
    # Add a "instrument_platform" column, currently not using this info,
    # just add this to be compatible with existing pipeline code
    if "instrument_platform" not in design.columns:
        design["instrument_platform"] = "ILLUMINA"

    # Add a "single_end" column
    design["single_end"] = design["read_2"].isna()

    # Add a "fasta" column, currently not using this info,
    # just add this to be compatible with existing pipeline code
    design["fasta"] = ""

    # Rename the read_1 and read_2 column to be consistent with existing code
    design = design.rename(columns={"read_1":"fastq_1", "read_2":"fastq_2"})

    # Extract group information to a different file
    metadata = design[["sample", "group"]]
    metadata.columns = ["sampleid", "group"]
    metadata.to_csv("group_metadata.csv", sep="\t", index=False)

    # Drop group column
    design = design.drop("group", axis=1)

    # Output
    design.to_csv(file_out, index=False)

def main(args=None):
    args = parse_args(args)
    check_samplesheet(args.FILE_IN, args.FILE_OUT)

if __name__ == "__main__":
    sys.exit(main())
