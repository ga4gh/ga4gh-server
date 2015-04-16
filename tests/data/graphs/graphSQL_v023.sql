-- Schema for submitting graphs to be served by GA4GH reference graph server
-- 
-- Graphs for the HGVM pilot bake-off can be loaded into the reference server using this
-- SQL schema. The file "graphData_v023.sql" shows an example dataset encoded in this format.
--
-- Contact maciek@soe.ucsc.edu for help with converting your graphs into this SQL format.
--
-- Fields not marked as NOT NULL are primarily optional metadata,
-- and their meaning is defined in corresponding GA4GH Avro records. 
--
-- FASTA URI table 
-- For current server version, just provide the filename in the URI field
--
CREATE TABLE FASTA (ID INTEGER PRIMARY KEY,
	fastaURI TEXT NOT NULL);
--
--
-- Sequences and joins - basis of graph topology
--
CREATE TABLE Sequence (ID INTEGER PRIMARY KEY,
	fastaID INTEGER NOT NULL REFERENCES FASTA(ID), -- the FASTA file that contains this sequence's bases.
	sequenceRecordName TEXT NOT NULL, -- access to the sequence bases in the FASTA file ONLY.
	md5checksum TEXT NOT NULL, -- checksum of the base sequence as found in the FASTA record.
	length INTEGER NOT NULL); -- length of the base sequence as found in the FASTA record.
--
--
CREATE TABLE GraphJoin (ID INTEGER PRIMARY KEY,
	--
	-- by convention, side1 < side2 in the lexicographic ordering defined by (sequenceID, position, forward).
	--
	side1SequenceID INTEGER NOT NULL REFERENCES Sequence(ID),
	side1Position INTEGER NOT NULL, -- 0 based indexing, counting from 5' end of sequence.
	side1StrandIsForward BOOLEAN NOT NULL, -- true if this side joins to 5' end of the base
	-- 
	side2SequenceID INTEGER NOT NULL REFERENCES Sequence(ID),
	side2Position INTEGER NOT NULL,
	side2StrandIsForward BOOLEAN NOT NULL);
--
--
-- References
--
CREATE TABLE Reference (ID INTEGER PRIMARY KEY, 
	name TEXT NOT NULL,
	updateTime DATE NOT NULL,
	sequenceID INTEGER NOT NULL REFERENCES Sequence(ID),
	start INTEGER, -- if null, reference starts at position 0 of the underlying sequence
	length INTEGER, -- if null, this is calculated as (sequence.lenght - start) 
	md5checksum TEXT, -- if null, assume sequence.md5checksum
	--
	-- the below metadata are defined as in the corresponding fields in the Avro Reference record.
	-- 
	isDerived BOOLEAN, 
	sourceDivergence REAL, 
	ncbiTaxonID INTEGER, 
	isPrimary BOOLEAN);
--
CREATE TABLE ReferenceAccession (ID INTEGER PRIMARY KEY,
	referenceID INTEGER NOT NULL REFERENCES Reference(ID),
	accessionID TEXT NOT NULL);
--
--
-- Reference sets
-- 
CREATE TABLE ReferenceSet (ID INTEGER PRIMARY KEY,
	ncbiTaxonID INT, -- may differ from ncbiTaxonID of contained Reference record
	description TEXT,
	assemblyID TEXT,
	isDerived BOOLEAN NOT NULL);
--
CREATE TABLE ReferenceSetAccession (ID INTEGER PRIMARY KEY,
	referenceSetID INTEGER NOT NULL REFERENCES ReferenceSet(ID),
	accessionID TEXT NOT NULL);
--
CREATE TABLE Reference_ReferenceSet_Join (referenceID INTEGER NOT NULL REFERENCES Reference(ID), 
	referenceSetID INTEGER NOT NULL REFERENCES ReferenceSet(ID),
	PRIMARY KEY(referenceID, referenceSetID));
--
CREATE TABLE GraphJoin_ReferenceSet_Join (graphJoinID INTEGER NOT NULL REFERENCES GraphJoin(ID),
	referenceSetID INTEGER NOT NULL REFERENCES ReferenceSet(ID),
	PRIMARY KEY(graphJoinID, referenceSetID));
--
--
-- Variant and call sets in the allelic world
--
CREATE TABLE VariantSet (ID INTEGER PRIMARY KEY,
	referenceSetID INTEGER NOT NULL REFERENCES ReferenceSet(ID),
	name TEXT);
--
CREATE TABLE CallSet (ID INTEGER PRIMARY KEY,
	name TEXT, -- can be null?
	sampleID TEXT);
--
CREATE TABLE VariantSet_CallSet_Join (variantSetID INTEGER NOT NULL REFERENCES VariantSet(ID), 
	callSetID INTEGER NOT NULL REFERENCES CallSet(ID),
	PRIMARY KEY(variantSetID, callSetID));
--
CREATE TABLE GraphJoin_VariantSet_Join (graphJoinID INTEGER NOT NULL REFERENCES GraphJoin(ID),
	variantSetID INTEGER NOT NULL REFERENCES VariantSet(ID),
	PRIMARY KEY(graphJoinID, variantSetID));
--
--
--
-- Allele and friends
--
CREATE TABLE Allele (ID INTEGER PRIMARY KEY, 
	variantSetID INTEGER REFERENCES VariantSet(ID), 
	name TEXT); -- Naming the allele is optional
--
CREATE TABLE AllelePathItem (alleleID INTEGER REFERENCES allele(ID), 
	pathItemIndex INTEGER NOT NULL, -- zero-based index of this pathItem within the entire path
	sequenceID INTEGER NOT NULL REFERENCES Sequence(ID), 
	start INTEGER NOT NULL,
	length INTEGER NOT NULL, 
	strandIsForward BOOLEAN NOT NULL,
	PRIMARY KEY(alleleID, pathItemIndex));
--
CREATE TABLE AlleleCall (alleleID INTEGER NOT NULL REFERENCES allele(ID), 
	callSetID INTEGER NOT NULL REFERENCES CallSet(ID),
	ploidy INTEGER NOT NULL,
	PRIMARY KEY(alleleID, callSetID));
--
--
