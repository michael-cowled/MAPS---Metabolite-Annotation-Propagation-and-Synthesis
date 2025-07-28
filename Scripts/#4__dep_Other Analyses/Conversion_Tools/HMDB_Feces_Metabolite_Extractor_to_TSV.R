require(xml2)
require(magrittr)
require(withr)
require(iterators)
require(itertools)
require(log4r)
require(methods)

InfoConsoleLv <- "DEBUG"

FileInputPath <- normalizePath("C:\\Temp\\feces_metabolites.xml" , mustWork = T)
FileOutputPath <- normalizePath("C:\\Temp\\feces_metabolites_target.TSV")

logObj <- logger(threshold = InfoConsoleLv ,
                 appenders = console_appender(layout = default_log_layout()))

openTag <- "<metabolite>"
closeTag <- "</metabolite>"

counter <- 0

Metabolite <- setRefClass(
  Class = "Metabolite" ,
  inheritPackage = T ,
  fields = list(
    accession = "character" ,
    hmdbStatus = "character" ,
    name = "character" ,
    smiles = "character" ,
    chemical_formula = "character" ,
    monisotopic_molecular_weight = "character" ,
    pubchem_compound_id = "character" ,
    XMLObj = "xml_document" ,
    normalConcentrationsCount = "numeric"
  ) ,
  
  methods = list(
    initialize = function(rawXMLText) {
      XMLObj <<- xml2::read_xml(rawXMLText)
      
      accession <<- extractAccession()
      
      hmdbStatus <<- extractStatus()
      
      name <<- extractName()
      smiles <<- extractSmiles()
      chemical_formula <<- extractChemicalFormula()
      monisotopic_molecular_weight <<- extractMonisotopicMolecularWeight()
      pubchem_compound_id <<- extractPubchemCompoundId()
      
      normalConcentrationsCount <<- extractNormalConcentration(countAll = T)
    } ,
    
    extractAccession = function(textJon = T , sep = ";") {
      primary <- xml_find_first(x = XMLObj , xpath = "./accession") %>%
        xml_text(trim = T)
      
      secondary <- xml_find_all(x = XMLObj , xpath = "./secondary_accessions/accession") %>%
        xml_text(trim = T)
      
      
      return(if (textJon) {
        paste(unique(c(primary , secondary)) , collapse = sep)
        
      } else{
        unique(c(primary , secondary))
      })
      
    } ,
    
    extractStatus = function() {
      xml_find_first(x = XMLObj , xpath = "./status") %>%
        xml_text(trim = T) %>%
        return()
    },
    
    extractName = function() {
      xml_find_first(x = XMLObj , xpath = "./name") %>%
        xml_text(trim = T) %>%
        return()
    } ,
    
    extractSmiles = function() {
      xml_find_first(x = XMLObj , xpath = "./smiles") %>%
        xml_text(trim = T) %>%
        return()
      
    } ,
    
    extractChemicalFormula = function() {
      xml_find_first(x = XMLObj , xpath = "./chemical_formula") %>%
        xml_text(trim = T) %>%
        return()
      
    } ,
    
    extractMonisotopicMolecularWeight = function() {
      xml_find_first(x = XMLObj , xpath = "./monisotopic_molecular_weight") %>%
        xml_text(trim = T) %>%
        return()
      
    } ,
    
    extractPubchemCompoundId = function() {
      xml_find_first(x = XMLObj , xpath = "./pubchem_compound_id") %>%
        xml_text(trim = T) %>%
        return()
      
    } ,
    
    extractNormalConcentration = function(countAll = T) {
      step1 <- xml2::xml_find_all(x = XMLObj , xpath = "./normal_concentrations/concentration[./biospecimen[text() = 'Feces']]")
      
      return(if (countAll) {
        length(step1)
      } else{
        step1
      })
      
    },
    
    formatOut = function(sep = "\t") {
      return(paste(
        c(
          accession,
          hmdbStatus,
          name,
          smiles,
          chemical_formula,
          monisotopic_molecular_weight,
          pubchem_compound_id,
          normalConcentrationsCount
        ) ,
        collapse = sep
      ))
      
    }
  )
  
)

Metabolite$accessors(
  c(
    "accession",
    "hmdbStatus",
    "name",
    "smiles",
    "chemical_formula",
    "monisotopic_molecular_weight",
    "pubchem_compound_id",
    "XMLObj",
    "normalConcentrationsCount"
  )
)

Metabolite$lock(
  c(
    "accession",
    "hmdbStatus",
    "name",
    "smiles",
    "chemical_formula",
    "monisotopic_molecular_weight",
    "pubchem_compound_id",
    "XMLObj",
    "normalConcentrationsCount"
  )
)


with_connection(con = list(
  fcon = file(FileInputPath , open = "rt") ,
  fout = file(FileOutputPath , open = "wt")
) , code = {
  writeLines(text =
               paste(
                 c(
                   "accession",
                   "name",
                   "smiles",
                   "chemical_formula",
                   "monisotopic_molecular_weight",
                   "pubchem_compound_id",
                   "normalConcentrationsCount"
                 ),
                 collapse = "\t"
               )  , con = fout)
  
  elemContainer <- NULL
  
  it <- ihasNext(ireadLines(con = fcon , n = 1))
  
  while (hasNext(it)) {
    linePtr <- nextElem(it)
    
    if (linePtr == openTag) {
      elemContainer <- vector(mode = "character")
      
      elemContainer %<>% append(linePtr)
      
    } else if (linePtr == closeTag) {
      elemContainer %<>% append(linePtr)
      
      counter %<>% add(1)
      
      metaboltePtr <- Metabolite$new(rawXMLText = paste(elemContainer , collapse = "\n"))
      
      outLinePtr <- metaboltePtr$formatOut()
      
      debug(logObj , outLinePtr)
      
      if (metaboltePtr$getNormalConcentrationsCount() > 0) {
        writeLines(text = outLinePtr , con = fout)
        
        for (nPtr in metaboltePtr$extractNormalConcentration(F)) {
          debug(logger = logObj , xml_text(
            xml_find_all(x = nPtr , xpath = ".//reference_text") ,
            trim = T
          ))
        }
      }
      
      # metaboltePtr$formatOut() %T>%
      #     debug( logger = logObj , counter ,  " --> " , . ) %>%
      #     writeLines(text = . , con = fout)
      
      if (!counter %% 100) {
        info(logger = logObj , "Processed " , counter , " metabolites")
      }
      
      elemContainer <- NULL
      
    } else {
      elemContainer %<>% append(linePtr)
    }
  }
})
