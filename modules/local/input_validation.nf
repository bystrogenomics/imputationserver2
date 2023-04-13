import groovy.json.JsonOutput

process INPUT_VALIDATION {

  publishDir params.output, mode: 'copy', pattern: '*.html'

  input:
    path(vcf_file)

  output:
    path("*.vcf.gz"), includeInputs: true, emit: validated_files
    path("*.html")
  script:

    config = [
        inputs: ['files'],
        params: [
            files: './',
            population: params.population,
            phasing: params.phasing,
            refpanel: params.refpanel.id,
            build: params.build,
            mode: params.mode
            //TODO: add missing params?
        ],
        data: [
            refpanel: params.refpanel
        ]
    ]

    """
    echo '${JsonOutput.toJson(config)}' > config.json

    java -cp /opt/imputationserver-utils/imputationserver-utils.jar \
      cloudgene.sdk.weblog.WebLogRunner \
      genepi.imputationserver.steps.InputValidation \
      config.json \
      01-input-validation.log

      ccat 01-input-validation.log --html > 01-input-validation.html

    """

}
