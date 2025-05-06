# Utility julia scripts

using CSV 
using DataFrames 
using Statistics
using SpecialFunctions # for cdf
using Distributions # for Normal distribution 

"""
    Function to calculate the p-value from a log file.
    It extracts the Total Observed scale h2 and its standard error from the log file,
    computes the z-score, and then calculates the p-value.
    
    Arguments:
    - log_file_path: String specifying the path to the log file
    
    Returns:
    - pval: Float64 representing the calculated p-value
    - z: Float64 representing the calculated z-score
"""
function calculate_heritability_p_value(log_file_path::String)

    log_text = read(log_file_path, String)
    pattern = r"Total Observed scale h2:\s+([0-9.]+)\s+\(([0-9.]+)\)"
    m = match(pattern, log_text)

    if m === nothing # OHHH nooo we couldn't find the pattern :(
        error("Pattern not found in log file.")
    end # NM my bad, we found the pattern

    h2 = parse(Float64, m.captures[1])
    se = parse(Float64, m.captures[2])

    z = h2 / se
    pval = 2 * (1 - cdf(Normal(0,1), abs(z)))

    println("z-score: ", z)
    println("p-value: %.10e\n", pval)

    return pval, z
end

"""
    Function to change the file extension of files in a specified folder.
    It reads each file, replaces spaces with tabs, and saves the updated content
    to a new file with a .tsv extension.
    
    Arguments:
    - datafolder: String specifying the path to the folder containing files
    - files: Vector of Strings specifying the names of files to be processed
    
    Returns:
    - None but prints "Done" when the operation is complete 
"""
function change_file_extension_for_empty_sep(files, datafolder) # see function change_file_extension in utility.jl
    for file in files
        input_path = joinpath(datafolder, file)

        output_filename = replace(file, ".txt" => ".tsv")
        output_path = joinpath(datafolder, output_filename)

        raw_text = read(input_path, String)
        cleaned_text = replace(raw_text, r"[ ]" => '\t')

        temp_path = joinpath(datafolder, "temp_cleaned_file.txt")
        
        open(temp_path, "w") do io
            write(io, cleaned_text)
        end

        df = CSV.File(temp_path; delim='\t', header=1) |> DataFrame
        CSV.write(output_path, df; delim='\t')
        rm(temp_path)
    end
end


function change_file_extension(file_list::Vector{String}) 
    for file in files
        input_path = joinpath(datafolder, file)

        output_filename = replace(file, ".txt" => ".tsv")
        output_path = joinpath(datafolder, output_filename)

        raw_text = read(input_path, String)
        cleaned_text = replace(raw_text, r"[ \t]+" => '\t')

        temp_path = joinpath(datafolder, "temp_cleaned_file.txt")
        
        open(temp_path, "w") do io
            write(io, cleaned_text)
        end

        df = CSV.File(temp_path; delim='\t', header=1) |> DataFrame
        CSV.write(output_path, df; delim='\t')
        rm(temp_path)
    end
end 


"""
    Function to fix column names in a TSV file by replacing spaces with underscores
    and removing parentheses. The updated DataFrame is saved to a specified output file.
    
    Arguments:
    - input_file: String specifying the path to the input TSV file
    - output_file: String specifying the path to save the updated TSV file
    
    Returns:
    - None but print "Done" when the operation is complete 
"""
function fix_column_names(input_file::String, output_file::String)
    df = CSV.read(input_file, DataFrame; delim='\t')
    new_names = [replace(replace(name, " " => "_"), "(" => "", ")" => "") for name in names(df)]
    rename!(df, Symbol.(new_names))
    CSV.write(output_file, df; delim='\t')
    println("Done")
end

"""
    Function to add a BETA column to a DataFrame based on the OR column.
    The BETA value is calculated as the natural logarithm of the OR value.
    The updated DataFrame is saved to a specified output file.
    
    Arguments:
    - df: DataFrame containing the data
    - output_file: String specifying the path to save the updated DataFrame
    
    Returns:
    - df: Updated DataFrame with the new BETA column
"""
function add_beta_and_save(df::DataFrame, output_file::String)
    
    if "OR" in names(df)
        df.BETA = log.(df.OR)
    end

    CSV.write(output_file, df; delim='\t')
    println("Updated DataFrame saved to: $output_file")
    
    return df
end

# The below function is NOT used in the pipeline but is useful for checking the headers of TSV files 
"""
    Function to compare headers of TSV files in a specified folder.
    It checks if the headers of all TSV files match the header of the first file.
    
    Arguments:
    - datafolder: String specifying the path to the folder containing TSV files
    
    Returns:
    - None but prints the comparison results
"""
function compare_headers(datafolder::String)

    tsv_files = filter(f -> endswith(f, ".tsv"), readdir(datafolder))

    for (i, file) in enumerate(tsv_files)

        reference_header = nothing
        file_path = joinpath(datafolder, file)
        println("\n--- File: $file ---")
        df = CSV.File(file_path; header=1, limit=0) |> DataFrame
        current_header = names(df)
        println("Headers: ", current_header)
        
        if i == 1
            reference_header = current_header
        else
            if current_header == reference_header
                println("Header matches reference.")
            else
                println("OHHH NOOO Header DOES NOT match reference!")
                println("Reference header: ", reference_header)
            end
        end
    end
end 

