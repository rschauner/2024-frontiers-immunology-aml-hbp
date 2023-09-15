

box::use(limma, methods[...])

.GetMetadata <- function(dge) {
    if (class(dge) == "ExpressionSet") {
        box::use(Biobase[phenoData])
        info <- phenoData(dge)@data
    }

    if (class(dge) == "DGEList") {
        info <- dge$samples
    }
    return(info)
}



AnnotationColors <- function() {
    list(
        age_group = c(Pediatric = "#011936", Adult = "#F9DC5C", Unknown = "#F4FFFD"),
        simple_disease = c(
            AML = "#426A5A", Early_Progenitors = "#DDAE7E", Whole_Blood = "#F2C57C", ALL = "#7FB685",
            Other = "#E9E6FF", Mixed = "#E9E6FF", MDS = "#E9E6FF", CML = "#E9E6FF"
        ),
        OGN = c(low = "#F6D8AE", OGA = "#004F2D", OGT = "#F4D35E", both = "#EF6F6C"),
        simple_disease2 = c(Adult = "#426A5A", Early_Progenitors = "#DDAE7E", Whole_Blood = "#F2C57C", Pediatric = "#7FB685")
    )
}

GetContrasts <- function(contrast, dge) {
    box::use(stats[...])
    info <- .GetMetadata(dge)

    hbp_design <- model.matrix(~ 0 + OGN, data = info)
    hbp_contrasts <- limma$makeContrasts(
        PvN = OGNboth - OGNlow,
        AvT = OGNOGA - OGNOGT,
        AvN = OGNOGA - OGNlow,
        TvN = OGNOGT - OGNlow,
        levels = c("OGNboth", "OGNlow", "OGNOGA", "OGNOGT")
    )

    disease_design <- model.matrix(~ 0 + simple_disease, data = info)
    disease_contrasts <- limma$makeContrasts(
        AMLvEP = simple_diseaseAML - simple_diseaseEarly_Progenitors,
        AMLvWB = simple_diseaseAML - simple_diseaseWhole_Blood,
        AMLvALL = simple_diseaseAML - simple_diseaseALL,
        levels = c(
            "simple_diseaseALL", "simple_diseaseAML", "simple_diseaseCML", "simple_diseaseMDS",
            "simple_diseaseMixed", "simple_diseaseOther",
            "simple_diseaseEarly_Progenitors", "simple_diseaseWhole_Blood"
        )
    )

    design_aml <- model.matrix(~ 0 + groups, data = info)
    contrasts_aml <- limma$makeContrasts(
        PvN = groupsboth_AML - groupslow_AML,
        AvT = groupsOGA_AML - groupsOGT_AML,
        AvN = groupsOGA_AML - groupslow_AML,
        TvN = groupsOGT_AML - groupslow_AML,
        levels = colnames(design_aml)
    )

    design_age <- model.matrix(~ 0 + age_with_disease, data = info)
    age_contrasts <- limma$makeContrasts(
        YvO = age_with_diseasePediatric_AML - age_with_diseaseAdult_AML,
        YvWB = age_with_diseasePediatric_AML - age_with_diseaseUnknown_Whole_Blood,
        OvWB = age_with_diseaseAdult_AML - age_with_diseaseUnknown_Whole_Blood,
        levels = colnames(design_age)
    )

    design_ratio <- model.matrix(~ 0 + OGT_OGA_ratio_bin, data = info)
    contrasts_ratio <- limma$makeContrasts(
        OGTvOGA = OGT_OGA_ratio_binOGT - OGT_OGA_ratio_binOGA,
        levels = colnames(design_ratio)
    )

    design_OGT <- model.matrix(~ 0 + OGT_bin, data = info)
    contrasts_OGT <- limma$makeContrasts(
        HvL = OGT_binOGT_high - OGT_binOGT_low,
        levels = colnames(design_OGT)
    )

    design_OGA <- model.matrix(~ 0 + OGA_bin, data = info)
    contrasts_OGA <- limma$makeContrasts(
        HvL = OGA_binOGA_high - OGA_binOGA_low,
        levels = colnames(design_OGA)
    )

    switch(
        contrast,
        HBP = list(hbp_design, hbp_contrasts),
        disease = list(disease_design, disease_contrasts),
        HBP_disease = list(design_aml, contrasts_aml),
        age = list(design_age, age_contrasts),
        ratio = list(design_ratio, contrasts_ratio),
        OGT = list(design_OGT, contrasts_OGT),
        OGA = list(design_OGA, contrasts_OGA)
    )
}

