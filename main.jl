# The julia script which runs on the whole pipeline 

using ArgParse 
using CSV
using DataFrames
using Conda

include("scripts/utility.jl")

function parse_commandline()
    s = ArgParseSettings() 

    "--rootfolder"
    help = "Path to the root folder"
    arg_type = String
    default = "-1"

    "--datafolder"
    help = "Path to the data folder"
    arg_type = String
    default = "data/sumstats"

    "--outfolder"
    help = "Path to the output folder"
    arg_type = String
    default = "data/results"

    "--ldsc_path"
    help = "Path to the ldsc folder"
    arg_type = String
    default = "ldsc-2.0.1" 

    "--reffolder"
    help = "Path to the reference folder"
    arg_type = String
    default = "data/reference"

    "--python_path"
    help = "Path to the python executable"
    arg_type = String
    default = "/Users/bingli/opt/miniconda3/envs/ldsc/bin/python"

    "--SUPERGNOVA_path"
    help = "Path to the SUPERGNOVA folder"
    arg_type = String
    default = "SUPERGNOVA"

    return parse_args(s) 

end
  
#-----------------------------------------------#       
#         Initialization
#-----------------------------------------------#
#--------------- Parse arguments ---------------# 
parsed_args = parse_commandline()
rootfolder = parsed_args["rootfolder"]
python_path = parsed_args["python_path"]
outfolder = parsed_args["outfolder"]
datafolder = parsed_args["datafolder"]
reffolder = parsed_args["reffolder"]
ldsc_path = parsed_args["ldsc_path"]
SUPERGNOVA_path = parsed_args["SUPERGNOVA_path"]

if rootfolder == "-1"
    println("Please provide the root folder path")
    rootfolder = pwd()
end

datafolder = joinpath(rootfolder, datafolder)
outfolder = joinpath(rootfolder, outfolder)
reffolder = joinpath(rootfolder, reffolder)
ldsc_path = joinpath(rootfolder, ldsc_path)
SUPERGNOVA_path = joinpath(rootfolder, SUPERGNOVA_path)
unused_data = joinpath(datafolder, "sumstats", "unused_data")
main_log_path = joinpath(rootfolder, "main.log")

mkpath(unused_data)

#--------------- Check data ---------------#
# compare_headers(datafolder) # okay seems that there are lost of mismatches for the headers 
# The below "check data" section is messy, but I will leave it here for now 

# change all .txt files to .tsv files 
files = ["adhd.txt", "income.txt", "baldness.txt", "cannabisuse.txt", "computeruse.txt"]

change_file_extension(files) # see function change_file_extension in utility.jl 

# income.txt is saved with [ ] instead of [\t]
# There must be much better way to code this 
# But this is faster for me ...   
files = ["income.txt"]
change_file_extension_for_empty_sep(files) # see function change_file_extension_for_empty_sep in utility.jl

chronotype_path = joinpath(datafolder, "chronotype.csv")
chronotype_outpath = joinpath(datafolder, "chronotype.tsv")
df = CSV.File(chronotype_path; delim=',', header=1) |> DataFrame 
CSV.write(chronotype_outpath, df; delim = "\t")
# Now all sumstates data have been converted to tsv files 

#--------------- heritability ---------------# 
# ---- Prepare data for running ldsc ---- # 

# ADHD data has no beta but only OR, need to add beta (lnOR) to the data 
# see function add_beta_and_save in utility.jl 
# Add log(OR) to ADHD and cannabis use 
adhd_path = joinpath(datafolder, "adhd.tsv")
df_adhd = CSV.File(adhd_path; delim='\t', header=1) |> DataFrame
modified_adhd = joinpath(datafolder, "adhd_modified.tsv")
add_beta_and_save(df_adhd, modified_adhd)

# move and rename folders
run(`mv $adhd_path $destination/`)
run(`mv $modified_adhd $/`)
# Later the function will be run on trait.tsv, so it's important to rename the file to trait.tsv
# and remove the original file to unused_data 

#-----------------------------------------------#

# Change column to make sure there is no space in the column names  
chronotype_path = joinpath(datafolder, "chronotype.tsv")
chronotype_outpath = joinpath(datafolder, "chronotype_cleaned.tsv")
# This function will remove the space in the column names
# see function fix_column_names in utility.jl 
# However, I found this function useless and a waste of time to write :-(
fix_column_names(chronotype_path, chronotype_outpath)
run(`mv $chronotype_path $unused_data`)
run(`mv $chronotype_outpath $chronotype_path`) 
# always make sure there are only one final trait.tsv file in the data folder

# double check the column names manually 
df = CSV.File(chronotype_path; delim='\t', header=1) |> DataFrame
println(names(df)) # ohhh nooo there is still a space in the column name 
rename!(df, Dict(
    Symbol("Meta\nln_OR") => :Meta_ln_OR
))
CSV.write(chronotype_path, df; delim='\t') # save the file again 

# computeruse data has NaN in the P column, impute with 0 
computeruse_path = joinpath(datafolder, "computeruse.tsv")
computeruse_outpath = joinpath(datafolder, "computeruse_cleaned.tsv")
df = CSV.File(computeruse_path; delim='\t', header=1) |> DataFrame

# This dataset has many empty columns and rows, need to remove them 
keep_cols = ["SNP", "CHR", "BP", "ALLELE1", "ALLELE0", "BETA", "P_BOLT_LMM_INF"] 
df = df[:, keep_cols] # only keep selected columns 
bad_rows = findall(row -> count(ismissing, row) > 0, eachrow(df)) # find all bad rows 
df = df[setdiff(1:nrow(df), bad_rows), :]

# check NaN values in the P column
df[!, "P_BOLT_LMM_INF"] = coalesce.(df[!, "P_BOLT_LMM_INF"], "NaN")
df[!, "P_BOLT_LMM_INF"] = replace.(df[!, "P_BOLT_LMM_INF"], r"0\.-" => "-")
df[!, "P_BOLT_LMM_INF"] = parse.(Float64, df[!, "P_BOLT_LMM_INF"]) 
if count(isnan, df[!, "P_BOLT_LMM_INF"]) != 0
    Error("Oh no! P column in computer use has NaN")# nothing -- good 
end

if count(isnan, df[!, "P_BOLT_LMM_INF"]) != 0 
    Error("Oh no! BETA column in computer use has NaN") # nothing -- good
end
CSV.write(computeruse_outpath, df; delim='\t') # save the file again
run(`mv $computeruse_path $unused_data`)
run(`mv $computeruse_outpath $computeruse_path`)

# Now all the data is cleaned and ready to be used 
# Quantitative traits dictionary 
quantitative_traits_info = Dict(

    "baldness" => Dict(
        "sample_size" => 205327,
        "columns" => Dict(
            "SNP" => "SNP",
            "CHR" => "CHR",
            "BP" => "BP",
            "ALLELE1" => "ALLELE1",
            "ALLELE0" => "ALLELE0",
            "BETA" => "BETA",
            "P" => "P_BOLT_LMM"
        )
    ),

    # "chronotype" => Dict(
    #     "sample_size" => 449734,
    #     "columns" => Dict(
    #         "SNP" => "Variant_ID",
    #         "CHR" => "CHR",
    #         "BP" => "POS",
    #         "ALLELE1" => "Allele1",
    #         "ALLELE0" => "Allele2",
    #         "BETA" => "Meta_ln_OR",
    #         "P" => "Meta_P_Corrected"
    #     )
    # ),

    "income" => Dict(
        "sample_size" => 286301,
        "columns" => Dict(
            "CHR" => "Chr",
            "SNP" => "SNP",
            "BP" => "BPos",
            "ALLELE0" => "Non_effect_Allele",
            "ALLELE1" => "Effect_Allele",
            "BETA" => "Beta",
            "P" => "P"
        )
    ),

    "computeruse" => Dict(
        "sample_size" => 408815,
        "columns" => Dict(
            "SNP" => "SNP",
            "CHR" => "CHR",
            "BP" => "BP",
            "ALLELE1" => "ALLELE1",
            "ALLELE0" => "ALLELE0",
            "BETA" => "BETA",
            "P" => "P_BOLT_LMM_INF"
        )
    )
)

# Binary traits dictionary with flipped columns
binary_traits_info = Dict(

    "ADHD" => Dict(
        "cases" => 38691,
        "controls" => 186843,
        "columns" => Dict(
            "SNP" => "SNP",
            "CHR" => "CHR",
            "BP" => "BP",
            "ALLELE1" => "A1",
            "ALLELE0" => "A2",
            "BETA" => "BETA",
            "P" => "P"
        )
    ),

    "cannabisuse" => Dict(
        "cases" => 14080,
        "controls" => 343726,
        "columns" => Dict(
            "SNP" => "SNP",
            "CHR" => "CHR",
            "BP" => "BP",
            "ALLELE1" => "A1",
            "ALLELE0" => "A2",
            "BETA" => "BETA",
            "P" => "P"
        )
    )
)


# run heritability analysis
function run_ldsc(datafolder::String, rootfolder::String, python_path::String, ref_path::String, output_dir, ldsc_path, trait_info::Dict, quantitative::Bool)

    mkpath(output_dir)  # make sure output directory exists

    for (trait, info) in trait_info

        println("Processing trait: $trait")

        if quantitative
            sample_size = info["sample_size"]
            bash_script_path = joinpath(rootfolder, "scripts/heritability_quantatitive.sh")
        else
            case_size = info["cases"] 
            control_size = info["controls"]
            bash_script_path = joinpath(rootfolder, "scripts/heritability_binary.sh")
        end

        required_fields = ["SNP", "ALLELE1", "ALLELE0", "P", "BETA"]
        columns = info["columns"]

        println("Columns map: ", columns)

        input_file = joinpath(datafolder, trait * ".tsv")
        output_prefix = joinpath(output_dir, trait * ".step1")

        # Check if all required fields exist
        for field in required_fields
            if !haskey(columns, field)
                error("Missing required key '$field' in columns for trait '$trait'!")
            end
        end

        SNP      = columns["SNP"]
        ALLELE1  = columns["ALLELE1"]
        ALLELE0  = columns["ALLELE0"]
        P        = columns["P"]
        BETA     = columns["BETA"]

        println("Mapped Columns: SNP=$SNP, ALLELE1=$ALLELE1, ALLELE0=$ALLELE0, P=$P, BETA=$BETA")

        println("Calling heritability.sh for trait: $trait")
        
        if quantitative
            run(`bash $bash_script_path \
                $datafolder $ldsc_path $sample_size $ref_path $output_dir \
                $SNP $ALLELE1 $ALLELE0 $P $BETA $trait $python_path`)
        else
            run(`bash $bash_script_path $datafolder $ldsc_path  \
            $case_size $control_size $ref_path $output_dir $SNP $ALLELE1 $ALLELE0  \
            $P $BETA $trait $python_path`)
        end 
    end
end

run_ldsc(datafolder, rootfolder, python_path, reffolder, outfolder, ldsc_path, quantitative_traits_info, true)
run_ldsc(datafolder, rootfolder, python_path, reffolder, outfolder, ldsc_path, binary_traits_info, false)

# ---- Calculate p-value for heritability ---- # 
baldness_step2_log = joinpath(outfolder, "baldness.step2.log")
open(main_log_path, "a") do io 
    p, z = calculate_heritability_p_value(baldness_step2_log)
    println(io, "Calculating p-value for heritability")
    println(io, "p-value for baldness: ", p)
    println(io, "z-score for baldness:", z)
end

#-----------------------------------------------# 
# ---- pairwise genetic correlation between baldness and 5 other traits ---- #
function run_genetic_correlation(rootfolder, python_path, reffolder, output_dir, ldsc_path, quantitative_traits_info, binary_traits_info)

    all_traits = collect(keys(quantitative_traits_info)) ∪ collect(keys(binary_traits_info))
    bash_script_path = joinpath(rootfolder, "scripts/pairwise_gencorrelation.sh")
    sumstats1 = joinpath(outfolder, "baldness.step1.sumstats.gz")

    for trait in all_traits

        println("Processing trait: $trait") 
    
        if trait != "baldness"
            sumstats2 = joinpath(outfolder, trait * ".step1.sumstats.gz")
            run(`bash $bash_script_path \
                $ldsc_path  $sumstats1 $sumstats2 \
                $reffolder $output_dir \
                $python_path "baldness" "$trait"`)
        end 
    end 
end 

run_genetic_correlation(rootfolder, python_path, reffolder, outfolder, ldsc_path, quantitative_traits_info, binary_traits_info)

# After running the bash script, we realized that chronotype have issues since the heritability is too low
# I couldn't find any better dataset with A1 and A2 columns and BETA althother
# The closest dataset is NUM of children saved in data/statsums but it doesn't have BETA 
# This cannot be simply transfered so let's do some 

numchildren_path = joinpath(datafolder, "numchildren.txt")
firstbirthmale_path = joinpath(datafolder, "firstbirthmale.txt") 
numchildren_outpath = joinpath(outfolder, "numchildren.step1.sumstats")
firstbirthmale_outpath = joinpath(outfolder, "firstbirthmale.step1.sumstats")

new_files = ["numchildren.txt", "firstbirthmale.txt"]
change_file_extension(new_files) # see function change_file_extension in utility.jl
# Okay a lot of missing values again, let's remove them 

function reformat_data(df_path::String, outpath::String, N::Int)
    df = CSV.File(df_path; delim='\t', header=1) |> DataFrame

    df_ldsc_ready = DataFrame(
    SNP = df.SNPID,
    A1 = df.A1,
    A2 = df.A2,
    Z = df.Zscore,
    N = fill(318463.000, nrow(df))
    )

    CSV.write(outpath, df_ldsc_ready, delim = "\t") # save the file again
    run(`gzip $outpath`) # compress the file
end


N_numchi = 318463 
N_firstbirthmale = 103736  

# reformat the data for numchildren and firstbirthmale 
reformat_data(numchildren_path, numchildren_outpath, N_numchi)
reformat_data(firstbirthmale_path, firstbirthmale_outpath, N_firstbirthmale) 

# Now we have the data ready for numchildren and firstbirthmale
# Let's re-run the genetic correlation analysis with these two traits between baldness
new_dict = Dict(
    "baldness" => Dict(
        "sample_size" => 205327
    ),
    "numchildren" => Dict(
        "sample_size" => 318463
    ),
    "firstbirthmale" => Dict(
        "sample_size" => 103736
    )
)

run_genetic_correlation(rootfolder, python_path, reffolder, outfolder, ldsc_path, new_dict, Dict())

#-----------Let's Do SUPERGNOVA -----------------# 
added_dict = Dict(
    "numchildren" => Dict(
        "sample_size" => 318463
    ),
    "firstbirthmale" => Dict(
        "sample_size" => 103736
    )
)

binary_size_dict = simplified_binary = Dict(
    trait => Dict("sample_size" => info["cases"] + info["controls"])
    for (trait, info) in binary_traits_info
)

quant_size_dic = Dict(
           trait => Dict("sample_size" => info["sample_size"])
           for (trait, info) in quantitative_traits_info
       )

sample_size_dict = merge(added_dict, binary_size_dict, quant_size_dic)


function run_supergnova_with_added_dict(sample_size_dict::Dict, rootfolder::String, outfolder::String, python_path::String, supergnova_path::String, reffolder::String)

    baldness_sumstats = joinpath(outfolder, "baldness.step1.sumstats.gz")
    N1 = sample_size_dict["baldness"]["sample_size"]
    bash_script_path = joinpath(rootfolder, "scripts/local_gencorrelation.sh")

    for (trait, info) in sample_size_dict
        if trait != "baldness"

            println("Running SuperGNOVA: baldness vs $trait")

            sumstats2 = joinpath(outfolder, trait * ".step1.sumstats.gz")
            N2 = info["sample_size"]

            bfile_path = joinpath(reffolder, "data", "bfiles", "eur_chr@_SNPmaf5")
            partition_path = joinpath(reffolder, "data", "partition", "eur_chr@.bed")

            output_path = joinpath(outfolder, "baldness_" * trait * ".supergnova.txt")
            
            # Let's just save the output in a new outfolder
            final_output = joinpath(outfolder, "local_genetic_correlation")
            mkpath(final_output) # make sure output directory exists

            run(`bash $bash_script_path \
                $baldness_sumstats $sumstats2 \
                $N1 $N2 \
                $bfile_path $partition_path \
                $final_output \
                $python_path $supergnova_path`)
        end
    end
end

# Remove any SNPs with missing values in the Z to run SUPERGNOVA 
files = [ "baldness.step1.sumstats.gz", "numchildren.step1.sumstats.gz", 
    "firstbirthmale.step1.sumstats.gz", "adhd.step1.sumstats.gz",
    "cannabisuse.step1.sumstats.gz", "computeruse.step1.sumstats.gz",
    "income.step1.sumstats.gz"]

for file in files
    path = joinpath(outfolder, file) 
    name = split(file, ".step1")[1]

    df = CSV.File(path; delim='\t', header=1) |> DataFrame
    outpath = joinpath(outfolder, name * "step1.sumstats")
    run(`mv $path $unused_data`) # move the file to unused_data folder

    # Remove rows with missing values in ALL columns
    filtered_df = dropmissing(df, ["N"])
    filtered_df = dropmissing(filtered_df, ["Z"])
    filtered_df = dropmissing(filtered_df, ["A1"])
    filtered_df = dropmissing(filtered_df, ["A2"])
    filtered_df = dropmissing(filtered_df, ["SNP"])

    CSV.write(outpath, filtered_df; delim='\t') 
    run(`gzip $outpath`) # compress the file
end 


run_supergnova_with_added_dict(sample_size_dict, rootfolder, outfolder, python_path, SUPERGNOVA_path, reffolder)


# ----------------- Run 



#----------------------- Plot making  ------------------------#
# ---- Build cell type specific plot ---- # 



