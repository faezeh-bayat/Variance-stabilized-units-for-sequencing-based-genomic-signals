input_directory="/absolute_path_of_the_users_input_files"
    # Options for "variance_stabilization_model" are "default" or "user_specified". If it is initialized to "default", VSS is going 
    # to use the default mean-varaince curve for stabilizing the variance of the user provided signals.
    # Otherwise, user need to provide two replicates for building the model.
    # variance_stabilization_model="default"
variance_stabilization_model="user_specified"
    # If user has specified the "user_specified" option for building the model, they should provide two replicate signal files.
    # bigWig, bedGraph and bam files are accepted for the provided signals.
replicate1_signals_for_training_the_model="rep1.bedGraph"
replicate2_signals_for_training_the_model="rep2.bedGraph"
chromosomes_to_build_the_model="chr21"

    # After the model is build, user may provide the signals they want to stabilize the variance. Same as building model procedure,
    # bigWig, bedGraph and bam files are accepted for the provided signals
signals_to_be_variance_stabilized="rep1.bedGraph"
chromosomes_to_be_stabilized="chr21"
Source_directory="/absolute_path_where_VSS_is_installed"

module load r
module load r/3.4.4
module load bedtools/2.27.1

cd "$Source_directory"
cd "Source"
Rscript VSS.R $input_directory \
                $variance_stabilization_model \
                $replicate1_signals_for_training_the_model \
                $replicate2_signals_for_training_the_model \
                $chromosomes_to_build_the_model \
                $signals_to_be_variance_stabilized \
                $chromosomes_to_be_stabilized \
                $Source_directory
                





