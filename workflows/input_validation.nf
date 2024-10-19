include { INPUT_VALIDATION_VCF } from '../modules/local/input_validation/input_validation_vcf'

// Handle wildcard paths and specific file paths
Channel
    .fromPath(params.files, checkIfExists: true)
    .flatMap { file -> 
        // Check if the path is a directory or a file
        def filePath = file.toString()
        if (file.isDirectory()) {
            // If it's a directory, find all .bgz and .gz files
            Channel.fromPath("${filePath}/*.{bgz,gz}")
        } else if (file.isFile()) {
            // If it's a file, just pass it along if it's a .vcf.gz, .bgz, or .gz file
            if (filePath.endsWith('.vcf.gz') || filePath.endsWith('.bgz') || filePath.endsWith('.gz')) {
                return [file]
            }
        }
    }
    .set { files }

workflow INPUT_VALIDATION {
    
    main:
    INPUT_VALIDATION_VCF(files.collect())

    emit:
    validated_files = INPUT_VALIDATION_VCF.out.validated_files.collect()
    validation_report = INPUT_VALIDATION_VCF.out.validation_report
}
