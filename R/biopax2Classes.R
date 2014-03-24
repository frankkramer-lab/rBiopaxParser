###############################################################################
#
# biopax2Classes.R: This file contains the class definitons, getters, setters
#					for the BioPax level 2 representation in R. This 
#					representation is based on S3-classes, so there is no real  
#					type-checking, inheritance, etc possible.
#					For more information on BioPax please visit biopax.org.
# author: Frank Kramer <dev@frankkramer.de>
#
# This is released under GPL-2.
# 
# Documentation was created using roxygen
#
###############################################################################

#' Class inheritance relationships in Biopax Level 2.
#' 
#' A data.frame listing all direct superclasses for every Biopax Level 2 class.
#' The variables are as follows:
#' 
#' \itemize{
#'   \item class. Name of the class
#'   \item superclass. Name of the superclass
#' }
#' 
#' @docType data
#' @keywords datasets
#' @name CLASS_INHERITANCE_BP2
#' @title CLASS_INHERITANCE_BP2
#' @usage CLASS_INHERITANCE_BP2
#' @format A data frame with 46 rows and 2 columns
#' @export
CLASS_INHERITANCE_BP2 = data.frame(
		matrix(ncol=2,byrow=T, dimnames=list(list(),list("class","superclass")),data= c(
"entity",								"",

"pathway",								"entity",

"interaction",							"entity",
"physicalInteraction",					"interaction",
"control",								"physicalInteraction",
"catalysis",							"control",
"modulation",							"control",
"conversion",							"physicalInteraction",
"complexAssembly",						"conversion",
"biochemicalReaction",					"conversion",
"transport",							"conversion",
"transportWithBiochemicalReaction",		"biochemicalReaction",
"transportWithBiochemicalReaction",		"transport",

"physicalEntity",						"entity",
"dna",									"physicalEntity",
"rna",									"physicalEntity",
"protein",								"physicalEntity",
"smallMolecule",						"physicalEntity",
"complex",								"physicalEntity",

"utilityClass",							"",

"chemicalStructure",					"utilityClass",
"deltaGprimeO",							"utilityClass",
"kPrime",								"utilityClass",
"confidence",							"utilityClass",
"evidence",								"utilityClass",
"experimentalForm",						"utilityClass",
"pathwayStep",							"utilityClass",
"sequenceFeature",						"utilityClass",
"sequenceLocation",						"utilityClass",

"sequenceInterval",						"sequenceLocation",
"sequenceSite",							"sequenceLocation",

"physicalEntityParticipant",			"utilityClass",

"sequenceParticipant",					"physicalEntityParticipant",
"dnaParticipant",						"physicalEntityParticipant",
"rnaParticipant",						"physicalEntityParticipant",
"proteinParticipant",					"physicalEntityParticipant",
"smallMoleculeParticipant",				"physicalEntityParticipant",
"complexParticipant",					"physicalEntityParticipant",

"externalReferenceUtilityClass",		"utilityClass",

"dataSource",							"externalReferenceUtilityClass",
"bioSource",							"externalReferenceUtilityClass",
"openControlledVocabulary",				"externalReferenceUtilityClass",
"xref",									"externalReferenceUtilityClass",

"unificationXref",						"xref",
"relationshipXref",						"xref",
"publicationXref",						"xref"
				)),	stringsAsFactors = FALSE
		)

#' Class properties in Biopax Level 2.
#' 
#' A data.frame listing all direct properties for every Biopax Level 2 class. 
#' Together with CLASS_INHERITANCE_BP2 this allows to list all properties, including the inherited ones, of every class.
#' 
#' The variables are as follows:
#' 
#' \itemize{
#'   \item class. Name of the class
#'   \item property. Name of the superclass
#'   \item property_type.Type of the property, value or reference
#'   \item cardinality. Maximum allowed cardinality of a property. Many properties may only be singular.
#' }
#' 
#' @docType data
#' @keywords datasets
#' @name CLASS_PROPERTIES_BP2
#' @title CLASS_PROPERTIES_BP2
#' @usage CLASS_PROPERTIES_BP2
#' @format A data frame with 106 rows and 4 columns
#' @export
CLASS_PROPERTIES_BP2 = data.frame(
		matrix(ncol=4,byrow=T, dimnames=list(list(),list("class","property","property_type","cardinality")),data= c(
"entity",								"NAME",						"http://www.w3.org/2001/XMLSchema#string",			"1",
"entity",								"SHORT-NAME",				"http://www.w3.org/2001/XMLSchema#string",			"1",
"entity",								"SYNONYMS",					"http://www.w3.org/2001/XMLSchema#string",			"*",
"entity",								"COMMENT",					"http://www.w3.org/2001/XMLSchema#string",			"*",
"entity",								"AVAILABILITY",				"http://www.w3.org/2001/XMLSchema#string",			"*",
"entity",								"DATA-SOURCE",				"dataSource",									"*",
"entity",								"XREF",						"xref",											"*",

"pathway",								"ORGANISM",					"bioSource",										"1",
"pathway",								"EVIDENCE",					"evidence",										"*",
"pathway",								"PATHWAY-COMPONENTS",		"interaction",									"*",
"pathway",								"PATHWAY-COMPONENTS",		"pathway",										"*",
"pathway",								"PATHWAY-COMPONENTS",		"pathwayStep",									"*",

"interaction",							"PARTICIPANTS",				"entity",										"*",
"interaction",							"PARTICIPANTS",				"physicalEntityParticipant",						"*",
"interaction",							"EVIDENCE",					"evidence",										"*",

"physicalInteraction",					"INTERACTION-TYPE",			"openControlledVocabulary",						"*",

"control",								"CONTROL-TYPE",				"http://www.w3.org/2001/XMLSchema#string",			"1",
"control",								"CONTROLLER",				"entity",										"*",
"control",								"CONTROLLER",				"physicalEntityParticipant",						"*",
"control",								"CONTROLLED",				"entity",										"*",
"control",								"CONTROLLED",				"physicalEntityParticipant",						"*",
"control",								"CONTROLLED",				"pathway",										"*",
"control",								"CONTROLLED",				"interaction",									"*",

"catalysis",							"DIRECTION",				"http://www.w3.org/2001/XMLSchema#string",			"1",
"catalysis",							"COFACTOR",					"entity",										"*",
"catalysis",							"COFACTOR",					"physicalEntityParticipant",						"*",

"conversion",							"SPONTANEUS",				"http://www.w3.org/2001/XMLSchema#string",			"1",
"conversion",							"LEFT",						"entity",										"*",
"conversion",							"LEFT",						"physicalEntityParticipant",						"*",
"conversion",							"RIGHT",					"entity",										"*",
"conversion",							"RIGHT",					"physicalEntityParticipant",						"*",

"biochemicalReaction",					"DELTA-H",					"http://www.w3.org/2001/XMLSchema#double",			"*",
"biochemicalReaction",					"DELTA-S",					"http://www.w3.org/2001/XMLSchema#double",			"*",
"biochemicalReaction",					"EC-NUMBER",				"http://www.w3.org/2001/XMLSchema#string",			"*",
"biochemicalReaction",					"DELTA-G",					"deltaGprimeO",									"*",
"biochemicalReaction",					"KEQ",						"kPrime",										"*",

"dna",									"SEQUENCE",					"http://www.w3.org/2001/XMLSchema#string",			"1",
"dna",									"ORGANISM",					"bioSource",										"1",
"rna",									"SEQUENCE",					"http://www.w3.org/2001/XMLSchema#string",			"1",
"rna",									"ORGANISM",					"bioSource",										"1",
"protein",								"SEQUENCE",					"http://www.w3.org/2001/XMLSchema#string",			"1",
"protein",								"ORGANISM",					"bioSource",										"1",

"complex",								"COMPONENTS",				"physicalEntityParticipant",						"*",
"complex",								"ORGANISM",					"bioSource",										"1",

"smallMolecule",						"MOLECULAR-WEIGHT",			"http://www.w3.org/2001/XMLSchema#double",			"1",
"smallMolecule",						"CHEMICAL-FORMULA",			"http://www.w3.org/2001/XMLSchema#string",			"1",
"smallMolecule",						"STRUCTURE",				"chemicalStructure",								"*",

"utilityClass",							"COMMENT",					"http://www.w3.org/2001/XMLSchema#string",			"*",

"chemicalStructure",					"STRUCTURE-DATA",			"http://www.w3.org/2001/XMLSchema#string",			"1",
"chemicalStructure",					"STRUCTURE-FORMAT",			"http://www.w3.org/2001/XMLSchema#string",			"1",

"confidence",							"XREF",						"publicationXref",								"*",
"confidence",							"CONFIDENCE-VALUE",			"http://www.w3.org/2001/XMLSchema#string",			"1",

"deltaGprimeO",							"DELTA-G-PRIME-O",			"http://www.w3.org/2001/XMLSchema#float",			"1",
"deltaGprimeO",							"IONIC-STRENGTH",			"http://www.w3.org/2001/XMLSchema#float",			"1",
"deltaGprimeO",							"PH",						"http://www.w3.org/2001/XMLSchema#float",			"1",
"deltaGprimeO",							"PMG",						"http://www.w3.org/2001/XMLSchema#float",			"1",
"deltaGprimeO",							"TEMPERATURE",				"http://www.w3.org/2001/XMLSchema#float",			"1",

"kPrime",								"IONIC-STRENGTH",			"http://www.w3.org/2001/XMLSchema#float",			"1",
"kPrime",								"PH",						"http://www.w3.org/2001/XMLSchema#float",			"1",
"kPrime",								"PMG",						"http://www.w3.org/2001/XMLSchema#float",			"1",
"kPrime",								"TEMPERATURE",				"http://www.w3.org/2001/XMLSchema#float",			"1",

"evidence",								"XREF",						"xref",											"*",
"evidence",								"CONFIDENCE",				"confidence",									"*",
"evidence",								"EVIDENCE-CODE",			"openControlledVocabulary",						"*",
"evidence",								"EXPERIMENTAL-FORM",		"experimentalForm",								"*",

"experimentalForm",						"EXPERIMENTAL-FORM-TYPE",	"openControlledVocabulary",						"*",
"experimentalForm",						"PARTICIPANT",				"physicalEntityParticipant",						"*",

"pathwayStep",							"NEXT-STEP",				"pathwayStep",									"*",
"pathwayStep",							"PATHWAY-COMPONENTS",		"interaction",									"*",
"pathwayStep",							"PATHWAY-COMPONENTS",		"pathway",										"*",
"pathwayStep",							"PATHWAY-COMPONENTS",		"pathwayStep",									"*",
"pathwayStep",							"STEP-INTERACTIONS",		"interaction",									"*",
"pathwayStep",							"STEP-INTERACTIONS",		"pathway",										"*",

"sequenceFeature",						"NAME",						"http://www.w3.org/2001/XMLSchema#string",			"1",
"sequenceFeature",						"SHORT-NAME",				"http://www.w3.org/2001/XMLSchema#string",			"1",
"sequenceFeature",						"SYNONYMS",					"http://www.w3.org/2001/XMLSchema#string",			"*",
"sequenceFeature",						"XREF",						"xref",											"*",
"sequenceFeature",						"FEATURE-TYPE",				"openControlledVocabulary",						"1",
"sequenceFeature",						"FEATURE-LOCATION",			"sequenceLocation",								"*",
"sequenceFeature",						"SEQUENCE-FEATURE-LIST",	"sequenceFeature",								"*",

"sequenceInterval",						"SEQUENCE-INTERVAL-BEGIN",	"sequenceSite",									"1",
"sequenceInterval",						"SEQUENCE-INTERVAL-END",	"sequenceSite",									"1",

"sequenceSite",							"POSITION-STATUS",			"http://www.w3.org/2001/XMLSchema#string",			"1",
"sequenceSite",							"SEQUENCE-POSITION",		"http://www.w3.org/2001/XMLSchema#integer",			"1",

"physicalEntityParticipant",			"STOICHIOMETRIC-COEFFICIENT","http://www.w3.org/2001/XMLSchema#double",			"1",
"physicalEntityParticipant",			"CELLULAR-LOCATION",		"openControlledVocabulary",						"1",
"physicalEntityParticipant",			"PHYSICAL-ENTITY",			"physicalEntity",								"1",

"sequenceParticipant",					"SEQUENCE-FEATURE-LIST",	"sequenceFeature",								"*",

"bioSource",							"NAME",						"http://www.w3.org/2001/XMLSchema#string",			"1",
"bioSource",							"CELLTYPE",					"openControlledVocabulary",						"1",
"bioSource",							"TISSUE",					"openControlledVocabulary",						"1",
"bioSource",							"TAXON-XREF",				"unificationXref",								"1",

"dataSource",							"NAME",						"http://www.w3.org/2001/XMLSchema#string",			"1",
"dataSource",							"XREF",						"xref",											"*",

"openControlledVocabulary",				"XREF",						"unificationXref",											"*",
"openControlledVocabulary",				"TERM",						"http://www.w3.org/2001/XMLSchema#string",			"*",

"xref",									"DB",						"http://www.w3.org/2001/XMLSchema#string",			"1",
"xref",									"DB-VERSION",				"http://www.w3.org/2001/XMLSchema#string",			"1",
"xref",									"ID",						"http://www.w3.org/2001/XMLSchema#string",			"1",
"xref",									"ID-VERSION",				"http://www.w3.org/2001/XMLSchema#string",			"1",

"relationshipXref",						"RELATIONSHIP-TYPE",		"http://www.w3.org/2001/XMLSchema#string",			"1",

"publicationXref",						"AUTHORS",					"http://www.w3.org/2001/XMLSchema#string",			"*",
"publicationXref",						"TITLE",					"http://www.w3.org/2001/XMLSchema#string",			"*",
"publicationXref",						"YEAR",						"http://www.w3.org/2001/XMLSchema#string",			"*",
"publicationXref",						"URL",						"http://www.w3.org/2001/XMLSchema#string",			"*",
"publicationXref",						"SOURCE",					"http://www.w3.org/2001/XMLSchema#string",			"*"

			)),	stringsAsFactors = FALSE		
		)

#' Class inheritance relationships in Biopax Level 3.
#' 
#' A data.frame listing all direct superclasses for every Biopax Level 3 class.
#' The variables are as follows:
#' 
#' \itemize{
#'   \item class. Name of the class
#'   \item superclass. Name of the superclass
#' }
#' 
#' NOT UPDATED TO BP3 yet!
#' 
#' @docType data
#' @keywords datasets
#' @name CLASS_INHERITANCE_BP3
#' @title CLASS_INHERITANCE_BP3
#' @usage CLASS_INHERITANCE_BP3
#' @format A data frame with 46 rows and 2 columns
#' @export
CLASS_INHERITANCE_BP3 = data.frame(
		matrix(ncol=2,byrow=T, dimnames=list(list(),list("class","superclass")),data= c(
						"Entity",								"",
						
						"Pathway",								"Entity",
						
						"Interaction",							"Entity",
						"PhysicalInteraction",					"Interaction",
						"Control",								"PhysicalInteraction",
						"Catalysis",							"Control",
						"Modulation",							"Control",
						"Conversion",							"PhysicalInteraction",
						"ComplexAssembly",						"Conversion",
						"BiochemicalReaction",					"Conversion",
						"Transport",							"Conversion",
						"TransportWithBiochemicalReaction",		"BiochemicalReaction",
						"TransportWithBiochemicalReaction",		"Transport",
						
						"PhysicalEntity",						"Entity",
						"Dna",									"PhysicalEntity",
						"Rna",									"PhysicalEntity",
						"Protein",								"PhysicalEntity",
						"SmallMolecule",						"PhysicalEntity",
						"Complex",								"PhysicalEntity",
						
						"UtilityClass",							"",
						
						"ChemicalStructure",					"UtilityClass",
						"DeltaG",					 			"UtilityClass",
						"kPrime",								"UtilityClass",
						"Evidence",								"UtilityClass",
						"ExperimentalForm",						"UtilityClass",
						"PathwayStep",							"UtilityClass",
						"EntityFeature",						"UtilityClass",
						"EntityRference",						"UtilityClass",
						"SequenceLocation",						"UtilityClass",
						
						"SequenceInterval",						"SequenceLocation",
						"SequenceSite",							"SequenceLocation",
						
						"PhysicalEntityParticipant",			"UtilityClass",
						
						"SequenceParticipant",					"PhysicalEntityParticipant",
						"DnaParticipant",						"PhysicalEntityParticipant",
						"RnaParticipant",						"PhysicalEntityParticipant",
						"ProteinParticipant",					"PhysicalEntityParticipant",
						"SmallMoleculeParticipant",				"PhysicalEntityParticipant",
						"ComplexParticipant",					"PhysicalEntityParticipant",
						
						"ExternalReferenceUtilityClass",		"UtilityClass",
						
						"BioSource",							"ExternalReferenceUtilityClass",
						"OpenControlledVocabulary",				"ExternalReferenceUtilityClass",
						"Xref",									"ExternalReferenceUtilityClass",
						
						"UnificationXref",						"Xref",
						"RelationshipXref",						"Xref",
						"PublicationXref",						"Xref"
				)),	stringsAsFactors = FALSE
)

#' Class properties in Biopax Level 3.
#' 
#' A data.frame listing all direct properties for every Biopax Level 3 class. 
#' Together with CLASS_INHERITANCE_BP3 this allows to list all properties, including the inherited ones, of every class.
#' 
#' The variables are as follows:
#' 
#' \itemize{
#'   \item class. Name of the class
#'   \item property. Name of the superclass
#'   \item property_type.Type of the property, value or reference
#'   \item cardinality. Maximum allowed cardinality of a property. Many properties may only be singular.
#' }
#' 
#' @docType data
#' @keywords datasets
#' @name CLASS_PROPERTIES_BP3
#' @title CLASS_PROPERTIES_BP3
#' @usage CLASS_PROPERTIES_BP3
#' @format A data frame with 106 rows and 4 columns
#' @export
CLASS_PROPERTIES_BP3 = data.frame(
		matrix(ncol=4,byrow=T, dimnames=list(list(),list("class","property","property_type","cardinality")),data= c(
						"Entity",								"name",						"http://www.w3.org/2001/XMLSchema#string",			"*",
						"Entity",								"displayName",				"http://www.w3.org/2001/XMLSchema#string",			"*",
						"Entity",								"standardName",				"http://www.w3.org/2001/XMLSchema#string",			"*",
						"Entity",								"comment",					"http://www.w3.org/2001/XMLSchema#string",			"*",
						"Entity",								"availability",				"http://www.w3.org/2001/XMLSchema#string",			"*",
						"Entity",								"dataSource",				"Provenance",									"*",
						"Entity",								"xref",						"Xref",											"*",
						
						"Pathway",								"organism",					"BioSource",										"*",
						"Pathway",								"evidence",					"Evidence",										"*",
						"Pathway",								"pathwayComponent",		"Interaction",									"*",
						"Pathway",								"pathwayComponent",		"Pathway",										"*",
						"Pathway",								"pathwayComponent",		"PathwayStep",									"*",
						
						"Interaction",							"participant",				"Entity",										"*",
						"Interaction",							"participant",				"physicalEntityParticipant",						"*",
						"Interaction",							"evidence",					"Evidence",										"*",
						"Interaction",							"interactionType",			"OpenControlledVocabulary",						"*",
						
						"Control",								"controlType",				"http://www.w3.org/2001/XMLSchema#string",			"*",
						"Control",								"controller",				"entity",										"*",
						"Control",								"controlled",				"entity",										"*",
						"Control",								"controlled",				"Pathway",										"*",
						"Control",								"controlled",				"Interaction",									"*",
						
						"Catalysis",							"direction",				"http://www.w3.org/2001/XMLSchema#string",			"*",
						"Catalysis",							"cofactor",					"entity",										"*",
						
						"Conversion",							"spontaneus",				"http://www.w3.org/2001/XMLSchema#string",			"*",
						"Conversion",							"left",						"entity",										"*",
						"Conversion",							"left",						"physicalEntityParticipant",						"*",
						"Conversion",							"right",					"entity",										"*",
						"Conversion",							"right",					"physicalEntityParticipant",						"*",
						
						"BiochemicalReaction",					"DELTA-H",					"http://www.w3.org/2001/XMLSchema#double",			"*",
						"BiochemicalReaction",					"DELTA-S",					"http://www.w3.org/2001/XMLSchema#double",			"*",
						"BiochemicalReaction",					"EC-NUMBER",				"http://www.w3.org/2001/XMLSchema#string",			"*",
						"BiochemicalReaction",					"DELTA-G",					"deltaGprimeO",									"*",
						"BiochemicalReaction",					"KEQ",						"kPrime",										"*",
						
						"Dna",									"sequence",					"http://www.w3.org/2001/XMLSchema#string",			"*",
						"Dna",									"organism",					"bioSource",										"*",
						"Rna",									"sequence",					"http://www.w3.org/2001/XMLSchema#string",			"*",
						"Rna",									"organism",					"bioSource",										"*",
						"Protein",								"sequence",					"http://www.w3.org/2001/XMLSchema#string",			"*",
						"Protein",								"organism",					"bioSource",										"*",
						
						"Complex",								"component",				"physicalEntityParticipant",						"*",
						"Complex",								"organism",					"bioSource",										"*",
						
						"SmallMolecule",						"molecularWeight",			"http://www.w3.org/2001/XMLSchema#double",			"*",
						"SmallMolecule",						"chemicalFormula",			"http://www.w3.org/2001/XMLSchema#string",			"*",
						"SmallMolecule",						"structure",				"chemicalStructure",								"*",
						
						"UtilityClass",							"comment",					"http://www.w3.org/2001/XMLSchema#string",			"*",
						
						"ChemicalStructure",					"structureData",			"http://www.w3.org/2001/XMLSchema#string",			"*",
						"ChemicalStructure",					"structureFormat",			"http://www.w3.org/2001/XMLSchema#string",			"*",
						
						"DeltaG",								"DELTA-G-PRIME-O",			"http://www.w3.org/2001/XMLSchema#float",			"*",
						"DeltaG",								"IONIC-STRENGTH",			"http://www.w3.org/2001/XMLSchema#float",			"*",
						"DeltaG",								"PH",						"http://www.w3.org/2001/XMLSchema#float",			"*",
						"DeltaG",								"PMG",						"http://www.w3.org/2001/XMLSchema#float",			"*",
						"DeltaG",								"temperature",				"http://www.w3.org/2001/XMLSchema#float",			"*",
						
						"kPrime",								"IONIC-STRENGTH",			"http://www.w3.org/2001/XMLSchema#float",			"*",
						"kPrime",								"PH",						"http://www.w3.org/2001/XMLSchema#float",			"*",
						"kPrime",								"PMG",						"http://www.w3.org/2001/XMLSchema#float",			"*",
						"kPrime",								"temperature",				"http://www.w3.org/2001/XMLSchema#float",			"*",
						
						"Evidence",								"xref",						"Xref",											"*",
						"Evidence",								"confidence",				"Confidence",									"*",
						"Evidence",								"evidenceCode",				"OpenControlledVocabulary",						"*",
						"Evidence",								"experimentalForm",			"ExperimentalForm",								"*",
						
						"ExperimentalForm",						"experimentalForm",			"OpenControlledVocabulary",						"*",
						"ExperimentalForm",						"participant",				"physicalEntityParticipant",						"*",
						
						"PathwayStep",							"nextStep",					"PathwayStep",									"*",
						"PathwayStep",							"pathwayComponent",			"Interaction",									"*",
						"PathwayStep",							"pathwayComponent",			"Pathway",										"*",
						"PathwayStep",							"pathwayComponent",			"PathwayStep",									"*",
						"PathwayStep",							"stepProcess",				"Interaction",									"*",
						"PathwayStep",							"stepProcess",				"Pathway",										"*",
						
						"EntityFeature",						"name",						"http://www.w3.org/2001/XMLSchema#string",			"*",
						"EntityFeature",						"SYNONYMS",					"http://www.w3.org/2001/XMLSchema#string",			"*",
						"EntityFeature",						"xref",						"Xref",											"*",
						"EntityFeature",						"FEATURE-TYPE",				"OpenControlledVocabulary",						"*",
						"EntityFeature",						"FEATURE-LOCATION",			"SequenceLocation",								"*",
						"EntityFeature",						"SEQUENCE-FEATURE-LIST",	"SequenceFeature",								"*",
						
						"SequenceInterval",						"SEQUENCE-INTERVAL-BEGIN",	"sequenceSite",									"*",
						"SequenceInterval",						"SEQUENCE-INTERVAL-END",	"sequenceSite",									"*",
						
						"SequenceSite",							"POSITION-STATUS",			"http://www.w3.org/2001/XMLSchema#string",			"*",
						"SequenceSite",							"SEQUENCE-POSITION",		"http://www.w3.org/2001/XMLSchema#integer",			"*",
						
						"BioSource",							"name",						"http://www.w3.org/2001/XMLSchema#string",			"*",
						"BioSource",							"CELLTYPE",					"OpenControlledVocabulary",						"*",
						"BioSource",							"TISSUE",					"OpenControlledVocabulary",						"*",
						"BioSource",							"TAXON-XREF",				"UnificationXref",								"*",
						
						"DataSource",							"name",						"http://www.w3.org/2001/XMLSchema#string",			"*",
						"DataSource",							"xref",						"Xref",											"*",
						
						"OpenControlledVocabulary",				"xref",						"UnificationXref",											"*",
						"OpenControlledVocabulary",				"TERM",						"http://www.w3.org/2001/XMLSchema#string",			"*",
						
						"Xref",									"DB",						"http://www.w3.org/2001/XMLSchema#string",			"*",
						"Xref",									"DB-VERSION",				"http://www.w3.org/2001/XMLSchema#string",			"*",
						"Xref",									"ID",						"http://www.w3.org/2001/XMLSchema#string",			"*",
						"Xref",									"ID-VERSION",				"http://www.w3.org/2001/XMLSchema#string",			"*",
						
						"RelationshipXref",						"RELATIONSHIP-TYPE",		"http://www.w3.org/2001/XMLSchema#string",			"*",
						
						"PublicationXref",						"AUTHORS",					"http://www.w3.org/2001/XMLSchema#string",			"*",
						"PublicationXref",						"TITLE",					"http://www.w3.org/2001/XMLSchema#string",			"*",
						"PublicationXref",						"YEAR",						"http://www.w3.org/2001/XMLSchema#string",			"*",
						"PublicationXref",						"URL",						"http://www.w3.org/2001/XMLSchema#string",			"*",
						"PublicationXref",						"SOURCE",					"http://www.w3.org/2001/XMLSchema#string",			"*"
				
				)),	stringsAsFactors = FALSE
)

		
#' This function returns the subclasses of the supplied biopax class.
#' 
#' This function returns the subclasses of the supplied biopax class.
#' 
#' @param classname A string containing a class name
#' @param biopaxlevel Numeric. Specifies the Biopax Level to use.
#' @return Returns character vector containing the subclasses of the supplied class
#' @author Frank Kramer
#' @export
#' @examples
#'  getSubClasses("control")
getSubClasses <- function(classname, biopaxlevel=3) {
	classname = stripns(classname)
	ret = list()
	for(x in 1:10) {
		if(biopaxlevel==2) {
			ret = CLASS_INHERITANCE_BP2$class[CLASS_INHERITANCE_BP2$superclass %in% c(classname,ret)]
		}
		if(biopaxlevel==3) {
			ret = CLASS_INHERITANCE_BP3$class[CLASS_INHERITANCE_BP3$superclass %in% c(classname,ret)]
		}
		
	}
	ret
}

#' This function returns the superclasses of the supplied biopax class.
#' 
#' This function returns the superclasses of the supplied biopax class.
#' 
#' @param classname A string containing a class name
#' @param biopaxlevel Numeric. Specifies the Biopax Level to use.
#' @return Returns character vector containing the superclasses of the supplied class
#' @author Frank Kramer
#' @export
#' @examples
#'  getSuperClasses("control")
getSuperClasses <- function(classname, biopaxlevel=3) {
	classname = stripns(classname)
	ret = list()
	for(x in 1:10) {
		if(biopaxlevel==2) {
			ret = CLASS_INHERITANCE_BP2$superclass[CLASS_INHERITANCE_BP2$class %in% c(classname,ret)]
		}
		if(biopaxlevel==3) {
			ret = CLASS_INHERITANCE_BP3$superclass[CLASS_INHERITANCE_BP3$class %in% c(classname,ret)]
		}
	}
	ret[ret != ""]
}

#' This function returns the properties of the supplied biopax class.
#' 
#' This function returns the properties of the supplied biopax class. It always considers inhertance.
#' Every class inhertis the properties of its super classes. A table listing all available properties and their cardinalities (for Biopax Level 2).
#' 
#' @param classname A string containing a class name
#' @param biopaxlevel Numeric. Specifies the Biopax Level to use.
#' @return Returns a data.frame containing the properties and cardinalities of the supplied class
#' @author Frank Kramer
#' @export
#' @examples
#'  getClassProperties("control")
getClassProperties <- function(classname, biopaxlevel=3) {
	classname = stripns(classname)
	if(biopaxlevel==2) {
		classes = c(classname,getSuperClasses(classname,biopaxlevel))
		return(CLASS_PROPERTIES_BP2[CLASS_PROPERTIES_BP2$class %in% classes,])
	}
	if(biopaxlevel==3) {
		classes = c(classname,getSuperClasses(classname,biopaxlevel))
		return(CLASS_PROPERTIES_BP3[CLASS_PROPERTIES_BP3$class %in% classes,])
	}
	return(NULL)
}

#
########### BIOPAX LEVEL 2 CORE CLASS DEFINITIONS
#
###########		ENTITY
##' 
##' @param instance_id	ID of the instance. Should be a a string, should be unique.
##' @param name 			
##' @param shortname 
##' @param synonyms 
##' @param comment 
##' @param datasource 
##' @param xref 
##' @param availability 
##' @returnType 
##' @return 
##' @author fkramer
##' @export
#newEntity <- function( instance_id, name, shortname=list(), synonyms=list(), comment=list(), datasource=list(),
#		xref=list(), availability=list() ){
#	ret = list(name=name, instance_id=instance_id,shortname=shortname,synonyms=synonyms,comment=comment,datasource=datasource,
#			xref=xref,availability=availability)
#	class(ret) <- c("biopax2_Entity",class(ret))
#	ret
#}
#
###########		PATHWAY
##' 
##' @param instance_id ID of the instance. Should be a a string, should be unique.
##' @param name 
##' @param shortname 
##' @param synonyms 
##' @param comment 
##' @param datasource 
##' @param xref 
##' @param availability 
##' @param organism 
##' @param evidence 
##' @param pathwaycomponents 
##' @returnType 
##' @return 
##' @author fkramer
##' @export
#newPathway <- function( instance_id, name, shortname=list(), synonyms=list(), comment=list(), datasource=list(),
#		xref=list(), availability=list(), organism=list(), evidence=list(), pathwaycomponents=list()){
#	ret = list(name=name, instance_id=instance_id,shortname=shortname,synonyms=synonyms,comment=comment,datasource=datasource,
#			xref=xref,availability=availability,organism=organism,evidence=evidence,pathwaycomponents=pathwaycomponents)
#	class(ret) <- c("biopax2_Pathway","biopax2_Entity",class(ret))
#	ret
#}
#
###########		PHYSICAL ENTITY
##' 
##' @param instance_id ID of the instance. Should be a a string, should be unique.
##' @param name 
##' @param shortname 
##' @param synonyms 
##' @param comment 
##' @param datasource 
##' @param xref 
##' @param availability 
##' @returnType 
##' @return 
##' @author fkramer
##' @export
#newPhysicalEntity <- function( instance_id, name, shortname=list(), synonyms=list(), comment=list(), datasource=list(),
#		xref=list(), availability=list() ){
#	ret = list(name=name, instance_id=instance_id,shortname=shortname,synonyms=synonyms,comment=comment,datasource=datasource,
#			xref=xref,availability=availability)
#	class(ret) <- c("biopax2_PhysicalEntity","biopax2_Entity",class(ret))
#	ret
#}
#
##' 
##' @param instance_id ID of the instance. Should be a a string, should be unique.
##' @param name 
##' @param shortname 
##' @param synonyms 
##' @param comment 
##' @param datasource 
##' @param xref 
##' @param availability 
##' @param sequence 
##' @param organism 
##' @returnType 
##' @return 
##' @author fkramer
##' @export
#newDNA <- function( instance_id, name, shortname=list(), synonyms=list(), comment=list(), datasource=list(),
#		xref=list(), availability=list(), sequence=list(), organism=list() ){
#	ret = list(name=name, instance_id=instance_id,shortname=shortname,synonyms=synonyms,comment=comment,datasource=datasource,
#			xref=xref,availability=availability, sequence=sequence, organism=organism)
#	class(ret) <- c("biopax2_DNA","biopax2_PhysicalEntity","biopax2_Entity",class(ret))
#	ret
#}
#
##' 
##' @param instance_id ID of the instance. Should be a a string, should be unique.
##' @param name 
##' @param shortname 
##' @param synonyms 
##' @param comment 
##' @param datasource 
##' @param xref 
##' @param availability 
##' @param sequence 
##' @param organism 
##' @returnType 
##' @return 
##' @author fkramer
##' @export
#newRNA <- function( instance_id, name, shortname=list(), synonyms=list(), comment=list(), datasource=list(),
#		xref=list(), availability=list(), sequence=list(), organism=list() ){
#	ret = list(name=name, instance_id=instance_id,shortname=shortname,synonyms=synonyms,comment=comment,datasource=datasource,
#			xref=xref,availability=availability, sequence=sequence, organism=organism)
#	class(ret) <- c("biopax2_RNA","biopax2_PhysicalEntity","biopax2_Entity",class(ret))
#	ret
#}
#
##' 
##' @param instance_id ID of the instance. Should be a a string, should be unique.
##' @param name 
##' @param shortname 
##' @param synonyms 
##' @param comment 
##' @param datasource 
##' @param xref 
##' @param availability 
##' @param sequence 
##' @param organism 
##' @returnType 
##' @return 
##' @author fkramer
##' @export
#newProtein <- function( instance_id, name, shortname=list(), synonyms=list(), comment=list(), datasource=list(),
#		xref=list(), availability=list(), sequence=list(), organism=list() ){
#	ret = list(name=name, instance_id=instance_id,shortname=shortname,synonyms=synonyms,comment=comment,datasource=datasource,
#			xref=xref,availability=availability, sequence=sequence, organism=organism)
#	class(ret) <- c("biopax2_Protein","biopax2_PhysicalEntity","biopax2_Entity",class(ret))
#	ret
#}
#
##' 
##' @param instance_id ID of the instance. Should be a a string, should be unique.
##' @param name 
##' @param shortname 
##' @param synonyms 
##' @param comment 
##' @param datasource 
##' @param xref 
##' @param availability 
##' @param components 
##' @param organism 
##' @returnType 
##' @return 
##' @author fkramer
##' @export
#newComplex <- function( instance_id, name, shortname=list(), synonyms=list(), comment=list(), datasource=list(),
#		xref=list(), availability=list(), components=list(), organism=list() ){
#	ret = list(name=name, instance_id=instance_id,shortname=shortname,synonyms=synonyms,comment=comment,datasource=datasource,
#			xref=xref,availability=availability, components=components, organism=organism)
#	class(ret) <- c("biopax2_Complex","biopax2_PhysicalEntity","biopax2_Entity",class(ret))
#	ret
#}
#
##' 
##' @param instance_id ID of the instance. Should be a a string, should be unique.
##' @param name 
##' @param shortname 
##' @param synonyms 
##' @param comment 
##' @param datasource 
##' @param xref 
##' @param availability 
##' @param structure 
##' @param chemicalformula 
##' @param molecularweight 
##' @returnType 
##' @return 
##' @author fkramer
##' @export
#newSmallMolecule <- function( instance_id, name, shortname=list(), synonyms=list(), comment=list(), datasource=list(),
#		xref=list(), availability=list(), structure=list(), chemicalformula=list(), molecularweight=list() ){
#	ret = list(name=name, instance_id=instance_id,shortname=shortname,synonyms=synonyms,comment=comment,datasource=datasource,
#			xref=xref,availability=availability,structure=structure,chemicalformula=chemicalformula,
#			molecularweight=molecularweight)
#	class(ret) <- c("biopax2_SmallMolecule","biopax2_PhysicalEntity","biopax2_Entity",class(ret))
#	ret
#}
#
###########		INTERACTION
##' 
##' @param instance_id ID of the instance. Should be a a string, should be unique.
##' @param name 
##' @param shortname 
##' @param synonyms 
##' @param comment 
##' @param datasource 
##' @param xref 
##' @param availability 
##' @param participants 
##' @param evidence 
##' @returnType 
##' @return 
##' @author fkramer
##' @export
#newInteraction <- function( instance_id, name, shortname=list(), synonyms=list(), comment=list(), datasource=list(),
#		xref=list(), availability=list(), participants=list(), evidence=list() ){
#	ret = list(name=name, instance_id=instance_id,shortname=shortname,synonyms=synonyms,comment=comment,datasource=datasource,
#			xref=xref,availability=availability,participants=participants,evidence=evidence)
#	class(ret) <- c("biopax2_Interaction","biopax2_Entity",class(ret))
#	ret
#}
#
##' 
##' @param instance_id ID of the instance. Should be a a string, should be unique.
##' @param name 
##' @param shortname 
##' @param synonyms 
##' @param comment 
##' @param datasource 
##' @param xref 
##' @param availability 
##' @param participants 
##' @param evidence 
##' @param interactiontype 
##' @returnType 
##' @return 
##' @author fkramer
##' @export
#newPhysicalInteraction <- function( instance_id, name, shortname=list(), synonyms=list(), comment=list(), datasource=list(),
#		xref=list(), availability=list(), participants=list(), evidence=list(), interactiontype=list() ){
#	ret = list(name=name, instance_id=instance_id,shortname=shortname,synonyms=synonyms,comment=comment,datasource=datasource,
#			xref=xref,availability=availability,participants=participants,evidence=evidence,
#			interactiontype=interactiontype)
#	class(ret) <- c("biopax2_PhysicalInteraction","biopax2_Interaction","biopax2_Entity",class(ret))
#	ret
#}
#
##' 
##' @param instance_id ID of the instance. Should be a a string, should be unique.
##' @param name 
##' @param shortname 
##' @param synonyms 
##' @param comment 
##' @param datasource 
##' @param xref 
##' @param availability 
##' @param participants 
##' @param evidence 
##' @param interactiontype 
##' @param controltype 
##' @param controller 
##' @param controlled 
##' @returnType 
##' @return 
##' @author fkramer
##' @export
#newControl <- function( instance_id, name, shortname=list(), synonyms=list(), comment=list(), datasource=list(),
#		xref=list(), availability=list(), participants=list(), evidence=list(), interactiontype=list(),
#		controltype=list(), controller=list(), controlled=list() ){
#	ret = list(name=name, instance_id=instance_id,shortname=shortname,synonyms=synonyms,comment=comment,datasource=datasource,
#			xref=xref,availability=availability,participants=participants,evidence=evidence,
#			interactiontype=interactiontype,controltype=controltype,controller=controller,controlled=controlled)
#	class(ret) <- c("biopax2_Control","biopax2_PhysicalInteraction","biopax2_Interaction","biopax2_Entity",class(ret))
#	ret
#}
#
##' 
##' @param instance_id ID of the instance. Should be a a string, should be unique.
##' @param name 
##' @param shortname 
##' @param synonyms 
##' @param comment 
##' @param datasource 
##' @param xref 
##' @param availability 
##' @param participants 
##' @param evidence 
##' @param interactiontype 
##' @param controltype 
##' @param controller 
##' @param controlled 
##' @param direction 
##' @param cofactor 
##' @returnType 
##' @return 
##' @author fkramer
##' @export
#newCatalysis <- function( instance_id, name, shortname=list(), synonyms=list(), comment=list(), datasource=list(),
#		xref=list(), availability=list(), participants=list(), evidence=list(), interactiontype=list(),
#		controltype=list(), controller=list(), controlled=list(), direction=list(), cofactor=list() ){
#	ret = list(name=name, instance_id=instance_id,shortname=shortname,synonyms=synonyms,comment=comment,datasource=datasource,
#			xref=xref,availability=availability,participants=participants,evidence=evidence,
#			interactiontype=interactiontype,controltype=controltype,controller=controller,controlled=controlled,
#			direction=direction,cofactor=cofactor)
#	class(ret) <- c("biopax2_Catalysis","biopax2_Control","biopax2_PhysicalInteraction","biopax2_Interaction","biopax2_Entity",class(ret))
#	ret
#}
#
##' 
##' @param instance_id ID of the instance. Should be a a string, should be unique.
##' @param name 
##' @param shortname 
##' @param synonyms 
##' @param comment 
##' @param datasource 
##' @param xref 
##' @param availability 
##' @param participants 
##' @param evidence 
##' @param interactiontype 
##' @param controltype 
##' @param controller 
##' @param controlled 
##' @param direction 
##' @param cofactor 
##' @returnType 
##' @return 
##' @author fkramer
##' @export
#newModulation <- function( instance_id, name, shortname=list(), synonyms=list(), comment=list(), datasource=list(),
#		xref=list(), availability=list(), participants=list(), evidence=list(), interactiontype=list(),
#		controltype=list(), controller=list(), controlled=list(), direction=list(), cofactor=list() ){
#	ret = list(name=name, instance_id=instance_id,shortname=shortname,synonyms=synonyms,comment=comment,datasource=datasource,
#			xref=xref,availability=availability,participants=participants,evidence=evidence,
#			interactiontype=interactiontype,controltype=controltype,controller=controller,controlled=controlled)
#	class(ret) <- c("biopax2_Modulation","biopax2_Control","biopax2_PhysicalInteraction","biopax2_Interaction","biopax2_Entity",class(ret))
#	ret
#}
#
##' 
##' @param instance_id ID of the instance. Should be a a string, should be unique.
##' @param name 
##' @param shortname 
##' @param synonyms 
##' @param comment 
##' @param datasource 
##' @param xref 
##' @param availability 
##' @param participants 
##' @param evidence 
##' @param interactiontype 
##' @param spontaneous 
##' @param left 
##' @param right 
##' @returnType 
##' @return 
##' @author fkramer
##' @export
#newConversion <- function( instance_id, name, shortname=list(), synonyms=list(), comment=list(), datasource=list(),
#		xref=list(), availability=list(), participants=list(), evidence=list(), interactiontype=list(),
#		spontaneous=list(), left=list(), right=list()){
#	ret = list(name=name, instance_id=instance_id,shortname=shortname,synonyms=synonyms,comment=comment,datasource=datasource,
#			xref=xref,availability=availability,participants=participants,evidence=evidence,
#			interactiontype=interactiontype,spontaneous=spontaneous,left=left,right=right)
#	class(ret) <- c("biopax2_Conversion","biopax2_PhysicalInteraction","biopax2_Interaction","biopax2_Entity",class(ret))
#	ret
#}
#
##' 
##' @param instance_id ID of the instance. Should be a a string, should be unique.
##' @param name 
##' @param shortname 
##' @param synonyms 
##' @param comment 
##' @param datasource 
##' @param xref 
##' @param availability 
##' @param participants 
##' @param evidence 
##' @param interactiontype 
##' @param spontaneous 
##' @param left 
##' @param right 
##' @returnType 
##' @return 
##' @author fkramer
##' @export
#newComplexAssembly <- function( instance_id, name, shortname=list(), synonyms=list(), comment=list(), datasource=list(),
#		xref=list(), availability=list(), participants=list(), evidence=list(), interactiontype=list(),
#		spontaneous=list(), left=list(), right=list()){
#	ret = list(name=name, instance_id=instance_id,shortname=shortname,synonyms=synonyms,comment=comment,datasource=datasource,
#			xref=xref,availability=availability,participants=participants,evidence=evidence,
#			interactiontype=interactiontype,spontaneous=spontaneous,left=left,right=right)
#	class(ret) <- c("biopax2_ComplexAssembly","biopax2_Conversion","biopax2_PhysicalInteraction","biopax2_Interaction","biopax2_Entity",class(ret))
#	ret
#}
#
##' 
##' @param instance_id ID of the instance. Should be a a string, should be unique.
##' @param name 
##' @param shortname 
##' @param synonyms 
##' @param comment 
##' @param datasource 
##' @param xref 
##' @param availability 
##' @param participants 
##' @param evidence 
##' @param interactiontype 
##' @param spontaneous 
##' @param left 
##' @param right 
##' @returnType 
##' @return 
##' @author fkramer
##' @export
#newTransport <- function( instance_id, name, shortname=list(), synonyms=list(), comment=list(), datasource=list(),
#		xref=list(), availability=list(), participants=list(), evidence=list(), interactiontype=list(),
#		spontaneous=list(), left=list(), right=list()){
#	ret = list(name=name, instance_id=instance_id,shortname=shortname,synonyms=synonyms,comment=comment,datasource=datasource,
#			xref=xref,availability=availability,participants=participants,evidence=evidence,
#			interactiontype=interactiontype,spontaneous=spontaneous,left=left,right=right)
#	class(ret) <- c("biopax2_Transport","biopax2_Conversion","biopax2_PhysicalInteraction","biopax2_Interaction","biopax2_Entity",class(ret))
#	ret
#}
#
##' 
##' @param instance_id ID of the instance. Should be a a string, should be unique.
##' @param name 
##' @param shortname 
##' @param synonyms 
##' @param comment 
##' @param datasource 
##' @param xref 
##' @param availability 
##' @param participants 
##' @param evidence 
##' @param interactiontype 
##' @param spontaneous 
##' @param left 
##' @param right 
##' @param deltag 
##' @param deltah 
##' @param deltas 
##' @param keq 
##' @param ecnumber 
##' @returnType 
##' @return 
##' @author fkramer
##' @export
#newBiochemicalReaction <- function( instance_id, name, shortname=list(), synonyms=list(), comment=list(), datasource=list(),
#		xref=list(), availability=list(), participants=list(), evidence=list(), interactiontype=list(),
#		spontaneous=list(), left=list(), right=list(), deltag=list(), deltah=list(), deltas=list(), keq=list(),
#		ecnumber=list()){
#	ret = list(name=name, instance_id=instance_id,shortname=shortname,synonyms=synonyms,comment=comment,datasource=datasource,
#			xref=xref,availability=availability,participants=participants,evidence=evidence,
#			interactiontype=interactiontype,spontaneous=spontaneous,left=left,right=right,
#			deltag=deltag,deltah=deltah,deltas=deltas,keq=keq,ecnumber=ecnumber)
#	class(ret) <- c("biopax2_BiochemicalReaction","biopax2_Conversion","biopax2_PhysicalInteraction","biopax2_Interaction","biopax2_Entity",class(ret))
#	ret
#}
#
##' 
##' @param instance_id ID of the instance. Should be a a string, should be unique.
##' @param name 
##' @param shortname 
##' @param synonyms 
##' @param comment 
##' @param datasource 
##' @param xref 
##' @param availability 
##' @param participants 
##' @param evidence 
##' @param interactiontype 
##' @param spontaneous 
##' @param left 
##' @param right 
##' @param deltag 
##' @param deltah 
##' @param deltas 
##' @param keq 
##' @param ecnumber 
##' @returnType 
##' @return 
##' @author fkramer
##' @export
#newTransportWithBiochemicalReaction <- function( instance_id, name, shortname=list(), synonyms=list(), comment=list(), datasource=list(),
#		xref=list(), availability=list(), participants=list(), evidence=list(), interactiontype=list(),
#		spontaneous=list(), left=list(), right=list(), deltag=list(), deltah=list(), deltas=list(), keq=list(),
#		ecnumber=list()){
#	ret = list(name=name, instance_id=instance_id,shortname=shortname,synonyms=synonyms,comment=comment,datasource=datasource,
#			xref=xref,availability=availability,participants=participants,evidence=evidence,
#			interactiontype=interactiontype,spontaneous=spontaneous,left=left,right=right,
#			deltag=deltag,deltah=deltah,deltas=deltas,keq=keq,ecnumber=ecnumber)
#	class(ret) <- c("biopax2_TransportWithBiochemicalReaction","biopax2_BiochemicalReaction","biopax2_Transport","biopax2_Conversion","biopax2_PhysicalInteraction","biopax2_Interaction","biopax2_Entity",class(ret))
#	ret
#}
#
############## END OF CORE CLASS DEFINITIONS
############## UTILITY CLASS DEFINITIONS
##utilityClass:
##chemicalStructure,confidence, evidence, externalReferenceUtilityClass, pathwayStep,
##physicalEntityParticipant, sequenceFeature, and sequenceLocation
#
###########		ENTITY
##' 
##' @param instance_id ID of the instance. Should be a a string, should be unique.
##' @param comment 
##' @returnType 
##' @return 
##' @author fkramer
##' @export
#newUtilityClass <- function( instance_id, comment=list() ){
#	ret = list(instance_id=instance_id,comment=comment)
#	class(ret) <- c("biopax2_UtilityClass",class(ret))
#	ret
#}
#
##' 
##' @param instance_id ID of the instance. Should be a a string, should be unique.
##' @param comment 
##' @param structureformat 
##' @param structuredata 
##' @returnType 
##' @return 
##' @author fkramer
##' @export
#newChemicalStructure <- function( instance_id, comment=list(), structureformat=list(), structuredata=list() ){
#	ret = list(instance_id=instance_id,comment=comment,structureformat=structureformat, structuredata=structuredata)
#	class(ret) <- c("biopax2_ChemicalStructure","biopax2_UtilityClass",class(ret))
#	ret
#}
#
##' 
##' @param instance_id ID of the instance. Should be a a string, should be unique.
##' @param comment 
##' @param confidencevalue 
##' @param xref 
##' @returnType 
##' @return 
##' @author fkramer
##' @export
#newConfidence <- function( instance_id, comment=list(), confidencevalue=list(), xref=list() ){
#	ret = list(instance_id=instance_id,comment=comment,confidencevalue=confidencevalue,xref=xref)
#	class(ret) <- c("biopax2_Confidence","biopax2_UtilityClass",class(ret))
#	ret
#}
#
##' 
##' @param instance_id ID of the instance. Should be a a string, should be unique.
##' @param comment 
##' @param xref 
##' @param experimentalform 
##' @param evidencecode 
##' @param confidence 
##' @returnType 
##' @return 
##' @author fkramer
##' @export
#newEvidence <- function( instance_id, comment=list(), xref=list(), experimentalform=list(), evidencecode=list(), confidence=list() ){
#	ret = list(instance_id=instance_id,comment=comment, xref=xref, experimentalform=experimentalform, evidencecode=evidencecode, confidence=confidence)
#	class(ret) <- c("biopax2_Evidence","biopax2_UtilityClass",class(ret))
#	ret
#}
#
##' 
##' @param instance_id ID of the instance. Should be a a string, should be unique.
##' @param comment 
##' @param experimentalformtype 
##' @param participant 
##' @returnType 
##' @return 
##' @author fkramer
##' @export
#newExperimentalForm <- function( instance_id, comment=list(), experimentalformtype=list(), participant=list() ){
#	ret = list(instance_id=instance_id,comment=comment, experimentalformtype=experimentalformtype, participant=participant)
#	class(ret) <- c("biopax2_ExperimentalForm","biopax2_UtilityClass",class(ret))
#	ret
#}
#
##' 
##' @param instance_id ID of the instance. Should be a a string, should be unique.
##' @param comment 
##' @param stepinteractions 
##' @param nextstep 
##' @returnType 
##' @return 
##' @author fkramer
##' @export
#newPathwayStep <- function( instance_id, comment=list(), stepinteractions=list(), nextstep=list() ){
#	ret = list(instance_id=instance_id,comment=comment, stepinteractions=stepinteractions, nextstep=nextstep)
#	class(ret) <- c("biopax2_PathwayStep","biopax2_UtilityClass",class(ret))
#	ret
#}
#
##' 
##' @param instance_id ID of the instance. Should be a a string, should be unique.
##' @param comment 
##' @param name 
##' @param featuretype 
##' @param shortname 
##' @param synonyms 
##' @param xref 
##' @param featurelocation 
##' @returnType 
##' @return 
##' @author fkramer
##' @export
#newSequenceFeature <- function( instance_id, comment=list(), name=list(), featuretype=list(), shortname=list(),
#								synonyms=list(), xref=list(), featurelocation=list() ){
#	ret = list(instance_id=instance_id,comment=comment, name=name, featuretype=featuretype, shortname=shortname,
#			synonyms=synonyms, xref=xref, featurelocation=featurelocation )
#	class(ret) <- c("biopax2_SequenceFeature","biopax2_UtilityClass",class(ret))
#	ret
#}
#
##' 
##' @param instance_id ID of the instance. Should be a a string, should be unique.
##' @param comment 
##' @returnType 
##' @return 
##' @author fkramer
##' @export
#newSequenceLocation <- function( instance_id, comment=list() ){
#	ret = list(instance_id=instance_id,comment=comment)
#	class(ret) <- c("biopax2_SequenceLocation","biopax2_UtilityClass",class(ret))
#	ret
#}
#
##' 
##' @param instance_id ID of the instance. Should be a a string, should be unique.
##' @param comment 
##' @param sequenceintervalbegin 
##' @param sequenceintervalend 
##' @returnType 
##' @return 
##' @author fkramer
##' @export
#newSequenceInterval <- function( instance_id, comment=list(), sequenceintervalbegin=list(), sequenceintervalend=list() ){
#	ret = list(instance_id=instance_id,comment=comment,sequenceintervalbegin=sequenceintervalbegin,sequenceintervalend=sequenceintervalend)
#	class(ret) <- c("biopax2_SequenceInterval","biopax2_SequenceLocation","biopax2_UtilityClass",class(ret))
#	ret
#}
#
##' 
##' @param instance_id ID of the instance. Should be a a string, should be unique.
##' @param comment 
##' @param positionstatus 
##' @param sequenceposition 
##' @returnType 
##' @return 
##' @author fkramer
##' @export
#newSequenceSite <- function( instance_id, comment=list(), positionstatus=list(), sequenceposition=list() ){
#	ret = list(instance_id=instance_id,comment=comment, positionstatus=positionstatus, sequenceposition=sequenceposition)
#	class(ret) <- c("biopax2_SequenceSite","biopax2_SequenceLocation","biopax2_UtilityClass",class(ret))
#	ret
#}
#
##' 
##' @param instance_id ID of the instance. Should be a a string, should be unique.
##' @param comment 
##' @param cellularlocation 
##' @param stoichiometriccoefficient 
##' @param physicalentity 
##' @returnType 
##' @return 
##' @author fkramer
##' @export
#newPhysicalEntityParticipant <- function( instance_id, comment=list(), cellularlocation=list(), stoichiometriccoefficient=list(), physicalentity=list() ){
#	ret = list(instance_id=instance_id,comment=comment, cellularlocation=cellularlocation, stoichiometriccoefficient=stoichiometriccoefficient, physicalentity=physicalentity)
#	class(ret) <- c("biopax2_PhysicalEntityParticipant","biopax2_UtilityClass",class(ret))
#	ret
#}
#
##' 
##' @param instance_id ID of the instance. Should be a a string, should be unique.
##' @param comment 
##' @param cellularlocation 
##' @param stoichiometriccoefficient 
##' @param physicalentity 
##' @param sequencefeaturelist 
##' @returnType 
##' @return 
##' @author fkramer
##' @export
#newSequenceParticipant <- function( instance_id, comment=list(), cellularlocation=list(), stoichiometriccoefficient=list(), physicalentity=list(), sequencefeaturelist=list() ){
#	ret = list(instance_id=instance_id,comment=comment, cellularlocation=cellularlocation, stoichiometriccoefficient=stoichiometriccoefficient, physicalentity=physicalentity, sequencefeaturelist=sequencefeaturelist)
#	class(ret) <- c("biopax2_SequenceParticipant","biopax2_PhysicalEntityParticipant","biopax2_UtilityClass",class(ret))
#	ret
#}
#
##' 
##' @param instance_id ID of the instance. Should be a a string, should be unique.
##' @param comment 
##' @param cellularlocation 
##' @param stoichiometriccoefficient 
##' @param physicalentity 
##' @returnType 
##' @return 
##' @author fkramer
##' @export
#newDNAParticipant <- function( instance_id, comment=list(), cellularlocation=list(), stoichiometriccoefficient=list(), physicalentity=list() ){
#	ret = list(instance_id=instance_id,comment=comment, cellularlocation=cellularlocation, stoichiometriccoefficient=stoichiometriccoefficient, physicalentity=physicalentity)
#	class(ret) <- c("biopax2_DNAParticipant","biopax2_PhysicalEntityParticipant","biopax2_UtilityClass",class(ret))
#	ret
#}
#
##' 
##' @param instance_id ID of the instance. Should be a a string, should be unique.
##' @param comment 
##' @param cellularlocation 
##' @param stoichiometriccoefficient 
##' @param physicalentity 
##' @returnType 
##' @return 
##' @author fkramer
##' @export
#newRNAParticipant <- function( instance_id, comment=list(), cellularlocation=list(), stoichiometriccoefficient=list(), physicalentity=list() ){
#	ret = list(instance_id=instance_id,comment=comment, cellularlocation=cellularlocation, stoichiometriccoefficient=stoichiometriccoefficient, physicalentity=physicalentity)
#	class(ret) <- c("biopax2_RNAParticipant","biopax2_PhysicalEntityParticipant","biopax2_UtilityClass",class(ret))
#	ret
#}
#
##' 
##' @param instance_id ID of the instance. Should be a a string, should be unique.
##' @param comment 
##' @param cellularlocation 
##' @param stoichiometriccoefficient 
##' @param physicalentity 
##' @returnType 
##' @return 
##' @author fkramer
##' @export
#newProteinEntityParticipant <- function( instance_id, comment=list(), cellularlocation=list(), stoichiometriccoefficient=list(), physicalentity=list() ){
#	ret = list(instance_id=instance_id,comment=comment, cellularlocation=cellularlocation, stoichiometriccoefficient=stoichiometriccoefficient, physicalentity=physicalentity)
#	class(ret) <- c("biopax2_ProteinParticipant","biopax2_PhysicalEntityParticipant","biopax2_UtilityClass",class(ret))
#	ret
#}
#
##' 
##' @param instance_id ID of the instance. Should be a a string, should be unique.
##' @param comment 
##' @param cellularlocation 
##' @param stoichiometriccoefficient 
##' @param physicalentity 
##' @returnType 
##' @return 
##' @author fkramer
##' @export
#newSmallMoleculeParticipant <- function( instance_id, comment=list(), cellularlocation=list(), stoichiometriccoefficient=list(), physicalentity=list() ){
#	ret = list(instance_id=instance_id,comment=comment, cellularlocation=cellularlocation, stoichiometriccoefficient=stoichiometriccoefficient, physicalentity=physicalentity)
#	class(ret) <- c("biopax2_SmallMoleculeParticipant","biopax2_PhysicalEntityParticipant","biopax2_UtilityClass",class(ret))
#	ret
#}
#
##' 
##' @param instance_id ID of the instance. Should be a a string, should be unique.
##' @param comment 
##' @param cellularlocation 
##' @param stoichiometriccoefficient 
##' @param physicalentity 
##' @returnType 
##' @return 
##' @author fkramer
##' @export
#newComplexParticipant <- function( instance_id, comment=list(), cellularlocation=list(), stoichiometriccoefficient=list(), physicalentity=list() ){
#	ret = list(instance_id=instance_id,comment=comment, cellularlocation=cellularlocation, stoichiometriccoefficient=stoichiometriccoefficient, physicalentity=physicalentity)
#	class(ret) <- c("biopax2_ComplexParticipant","biopax2_PhysicalEntityParticipant","biopax2_UtilityClass",class(ret))
#	ret
#}
#
##' 
##' @param instance_id ID of the instance. Should be a a string, should be unique.
##' @param comment 
##' @returnType 
##' @return 
##' @author fkramer
##' @export
#newExternalReferenceUtilityClass <- function( instance_id, comment=list() ){
#	ret = list(instance_id=instance_id,comment=comment)
#	class(ret) <- c("biopax2_ExternalReferenceUtilityClass","biopax2_UtilityClass",class(ret))
#	ret
#}
#
##' 
##' @param instance_id ID of the instance. Should be a a string, should be unique.
##' @param comment 
##' @param name 
##' @param xref 
##' @returnType 
##' @return 
##' @author fkramer
##' @export
#newDataSource <- function( instance_id, comment=list(), name=list(), xref=list() ){
#	ret = list(instance_id=instance_id,comment=comment, name=name, xref=xref)
#	class(ret) <- c("biopax2_DataSource","biopax2_ExternalReferenceUtilityClass","biopax2_UtilityClass",class(ret))
#	ret
#}
#
##' 
##' @param instance_id ID of the instance. Should be a a string, should be unique.
##' @param comment 
##' @param name 
##' @param taxonxref 
##' @param tissue 
##' @param celltype 
##' @returnType 
##' @return 
##' @author fkramer
##' @export
#newBioSource <- function( instance_id, comment=list(), name=list(), taxonxref=list(), tissue=list(), celltype=list() ){
#	ret = list(instance_id=instance_id,comment=comment, name=name, taxonxref=taxonxref, tissue=tissue, celltype=celltype)
#	class(ret) <- c("biopax2_BioSource","biopax2_ExternalReferenceUtilityClass","biopax2_UtilityClass",class(ret))
#	ret
#}
#
##' 
##' @param instance_id ID of the instance. Should be a a string, should be unique.
##' @param comment 
##' @param term 
##' @param xref 
##' @returnType 
##' @return 
##' @author fkramer
##' @export
#newOpenControlledVocabulary <- function( instance_id, comment=list(), term=list(), xref=list() ){
#	ret = list(instance_id=instance_id,comment=comment, term=term, xref=xref)
#	class(ret) <- c("biopax2_OpenControlledVocabulary","biopax2_ExternalReferenceUtilityClass","biopax2_UtilityClass",class(ret))
#	ret
#}
#
##' 
##' @param instance_id ID of the instance. Should be a a string, should be unique.
##' @param comment 
##' @param db 
##' @param dbversion 
##' @param id 
##' @param idversion 
##' @returnType 
##' @return 
##' @author fkramer
##' @export
#newXref <- function( instance_id, comment=list(), db=list(), dbversion=list(), id=list(), idversion=list() ){
#	ret = list(instance_id=instance_id,comment=comment, db=db, dbversion=dbversion, id=id, idversion=idversion)
#	class(ret) <- c("biopax2_Xref","biopax2_ExternalReferenceUtilityClass","biopax2_UtilityClass",class(ret))
#	ret
#}
#
##' 
##' @param instance_id ID of the instance. Should be a a string, should be unique.
##' @param comment 
##' @param db 
##' @param dbversion 
##' @param id 
##' @param idversion 
##' @returnType 
##' @return 
##' @author fkramer
##' @export
#newUnificationXref <- function( instance_id, comment=list(), db=list(), dbversion=list(), id=list(), idversion=list() ){
#	ret = list(instance_id=instance_id,comment=comment, db=db, dbversion=dbversion, id=id, idversion=idversion)
#	class(ret) <- c("biopax2_UnificationXref","biopax2_Xref","biopax2_ExternalReferenceUtilityClass","biopax2_UtilityClass",class(ret))
#	ret
#}
#
##' 
##' @param instance_id ID of the instance. Should be a a string, should be unique.
##' @param comment 
##' @param db 
##' @param dbversion 
##' @param id 
##' @param idversion 
##' @param relationshiptype 
##' @returnType 
##' @return 
##' @author fkramer
##' @export
#newRelationshipXref <- function( instance_id, comment=list(), db=list(), dbversion=list(), id=list(), idversion=list(), relationshiptype=list() ){
#	ret = list(instance_id=instance_id,comment=comment, db=db, dbversion=dbversion, id=id, idversion=idversion, relationshiptype=relationshiptype)
#	class(ret) <- c("biopax2_RelationshipXref","biopax2_Xref","biopax2_ExternalReferenceUtilityClass","biopax2_UtilityClass",class(ret))
#	ret
#}
#
##' 
##' @param instance_id ID of the instance. Should be a a string, should be unique.
##' @param comment 
##' @param db 
##' @param dbversion 
##' @param id 
##' @param idversion 
##' @param title 
##' @param year 
##' @param authors 
##' @param url 
##' @param source 
##' @returnType 
##' @return 
##' @author fkramer
##' @export
#newPublicationXref <- function( instance_id, comment=list(), db=list(), dbversion=list(), id=list(), idversion=list(), title=list(), year=list(), authors=list(), url=list(), source=list() ){
#	ret = list(instance_id=instance_id,comment=comment, db=db, dbversion=dbversion, id=id, idversion=idversion, title=title, year=year, authors=authors, url=url, source=source)
#	class(ret) <- c("biopax2_PublicationXref","biopax2_Xref","biopax2_ExternalReferenceUtilityClass","biopax2_UtilityClass",class(ret))
#	ret
#}
#
############## END OF CLASS DEFINITIONS
