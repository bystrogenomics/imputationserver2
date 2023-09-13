if (params.refpanel_yaml){
    params.refpanel = RefPanelUtil.loadFromFile(params.refpanel_yaml)
}

include { NO_PHASING               } from '../modules/local/no_phasing'
include { PHASING_EAGLE            } from '../modules/local/phasing_eagle'
include { PHASING_BEAGLE           } from '../modules/local/phasing_beagle'

workflow PHASING_WF {

    take: 
        chunks_vcf 
        chunks_csv
    main:

    chunks_vcf
        .flatten()
        .map { it -> tuple(file(it).baseName, it) }
        .set{ chunks_vcf_index }

    chunks_csv
        .flatten()
        .splitCsv(header:false, sep:'\t')
        .map{ 
            row-> tuple(file(row[4]).baseName, row[0], row[1], row[2], row[3], row[4], row[5], row[6])
        }
        .set { chunks_csv_index }

    chunks_csv_index
        .combine(chunks_vcf_index, by: 0)
        .map{
            row-> tuple(row[1], row[2], row[3], row[4], file(row[8]), row[6], row[7])
        }
        .set { metafiles_ch }


    if ("${params.refpanel.refEagle}" != null) {

        autosomes_eagle_ch =  Channel.from ( 1..22)
        .map { it -> tuple(it.toString(), file("$params.refpanel.refEagle".replaceAll('\\$chr', it.toString())),file("$params.refpanel.refEagle".replaceAll('\\$chr', it.toString())+'.csi')) }

        non_autosomes_eagle_ch =  Channel.from ( 'X.nonPAR', 'X.PAR1', 'X.PAR2', 'MT')
        .map { it -> tuple(it.toString(), file("$params.refpanel.refEagle".replaceAll('\\$chr', it.toString())),file("$params.refpanel.refEagle".replaceAll('\\$chr', it.toString())+'.csi')) }

        eagle_bcf_ch = autosomes_eagle_ch.concat(non_autosomes_eagle_ch)

    }

    if ("${params.refpanel.refBeagle}" != null) {

        autosomes_beagle_ch = Channel.from ( 1..22 )
        .map { it -> tuple(it.toString(), file("$params.refpanel.refBeagle".replaceAll('\\$chr', it.toString()))) }

        non_autosomes_beagle_ch = Channel.from ( 'X.nonPAR', 'X.PAR1', 'X.PAR2', 'MT')
        .map { it -> tuple(it.toString(), file("$params.refpanel.refBeagle".replaceAll('\\$chr', it.toString()))) }

        beagle_bcf_ch = autosomes_beagle_ch.concat(non_autosomes_beagle_ch)

        autosomes_beagle_map_ch = Channel.from ( 1..22 )
        .map { it -> tuple(it.toString(), file("$params.refpanel.mapBeagle".replaceAll('\\$chr', it.toString()))) }

        non_autosomes_beagle_map_ch = Channel.from (  'X.nonPAR', 'X.PAR1', 'X.PAR2', 'MT' )
        .map { it -> tuple(it.toString(), file("$params.refpanel.mapBeagle".replaceAll('\\$chr', it.toString()))) }

        beagle_map_ch = autosomes_beagle_map_ch.concat(non_autosomes_beagle_map_ch)
    }

    autosomes_m3vcf_ch = Channel.from ( 1..22 )
        .map { it -> tuple(it.toString(), file("$params.refpanel.hdfs".replaceAll('\\$chr', it.toString()))) }

    non_autosomes_m3vcf_ch = Channel.from ( 'X.nonPAR', 'X.PAR1', 'X.PAR2', 'MT')
        .map { it -> tuple(it.toString(), file("$params.refpanel.hdfs".replaceAll('\\$chr', it.toString()))) }

    minimac_m3vcf_ch = autosomes_m3vcf_ch.concat(non_autosomes_m3vcf_ch)

    //TODO: read from Dockerfile
    // check for '' required for testPipelineWithPhasedAndEmptyPhasing. Test case could be deleted since phasing is never '' anymore
    if ("${params.phasing}" == 'eagle' || "${params.phasing}" == '') {
    phasing_method = params.eagle_version
    }
    else if ("${params.phasing}" == 'beagle') {
    phasing_method = params.beagle_version
    }
    else if ("${params.phasing}" == 'no_phasing') {
    phasing_method = "n/a"
    }

    map_eagle   = file(params.refpanel.mapEagle, checkIfExists: false)
    //map_beagle  = file(params.refpanel.mapBeagle, checkIfExists: false)


    // check for '' required for testPipelineWithPhasedAndEmptyPhasing. Test case could be deleted since phasing is never '' anymore
    if ("${params.phasing}" == 'eagle'  || "${params.phasing}" == '') {

     eagle_bcf_metafiles_ch =  eagle_bcf_ch.combine(metafiles_ch, by: 0)

     PHASING_EAGLE ( eagle_bcf_metafiles_ch, map_eagle, phasing_method )

     phased_m3vcf_ch = PHASING_EAGLE.out.eagle_phased_ch.combine(minimac_m3vcf_ch, by: 0)

    }

    if ("${params.phasing}" == 'beagle') {

     beagle_bcf_metafiles_ch = beagle_bcf_ch.combine(metafiles_ch, by: 0)

     //combine with map since also split by chromsome
     beagle_bcf_metafiles_map_ch = beagle_bcf_metafiles_ch.combine(beagle_map_ch, by: 0)

     PHASING_BEAGLE ( beagle_bcf_metafiles_map_ch, phasing_method )

     phased_m3vcf_ch = PHASING_BEAGLE.out.beagle_phased_ch.combine(minimac_m3vcf_ch, by: 0)

    }

    if ("${params.phasing}" == 'no_phasing') {

     NO_PHASING (metafiles_ch)

     phased_m3vcf_ch = NO_PHASING.out.skipped_phasing_ch.combine(minimac_m3vcf_ch, by: 0)

    }

    emit: phased_m3vcf_ch

}


workflow.onComplete {
    println "Pipeline completed at: $workflow.complete"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}


