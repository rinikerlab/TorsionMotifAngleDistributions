--
-- PostgreSQL database dump
--

-- Dumped from database version 16.9 (Ubuntu 16.9-1.pgdg22.04+1)
-- Dumped by pg_dump version 16.9 (Ubuntu 16.9-1.pgdg22.04+1)

SET statement_timeout = 0;
SET lock_timeout = 0;
SET idle_in_transaction_session_timeout = 0;
SET client_encoding = 'UTF8';
SET standard_conforming_strings = on;
SELECT pg_catalog.set_config('search_path', '', false);
SET check_function_bodies = false;
SET xmloption = content;
SET client_min_messages = warning;
SET row_security = off;

--
-- Name: bookkeeping; Type: SCHEMA; Schema: -; Owner: USERNAMEPLACEHOLDER
--

CREATE SCHEMA bookkeeping;


ALTER SCHEMA bookkeeping OWNER TO USERNAMEPLACEHOLDER;

--
-- Name: comparisons; Type: SCHEMA; Schema: -; Owner: USERNAMEPLACEHOLDER
--

CREATE SCHEMA comparisons;


ALTER SCHEMA comparisons OWNER TO USERNAMEPLACEHOLDER;

--
-- Name: csd202403; Type: SCHEMA; Schema: -; Owner: USERNAMEPLACEHOLDER
--

CREATE SCHEMA csd202403;


ALTER SCHEMA csd202403 OWNER TO USERNAMEPLACEHOLDER;

--
-- Name: csd202403modified; Type: SCHEMA; Schema: -; Owner: USERNAMEPLACEHOLDER
--

CREATE SCHEMA csd202403modified;


ALTER SCHEMA csd202403modified OWNER TO USERNAMEPLACEHOLDER;

--
-- Name: csd202403modifiedrdk; Type: SCHEMA; Schema: -; Owner: USERNAMEPLACEHOLDER
--

CREATE SCHEMA csd202403modifiedrdk;


ALTER SCHEMA csd202403modifiedrdk OWNER TO USERNAMEPLACEHOLDER;

--
-- Name: csd202403rdk; Type: SCHEMA; Schema: -; Owner: USERNAMEPLACEHOLDER
--

CREATE SCHEMA csd202403rdk;


ALTER SCHEMA csd202403rdk OWNER TO USERNAMEPLACEHOLDER;

--
-- Name: extractions; Type: SCHEMA; Schema: -; Owner: USERNAMEPLACEHOLDER
--

CREATE SCHEMA extractions;


ALTER SCHEMA extractions OWNER TO USERNAMEPLACEHOLDER;

--
-- Name: fits; Type: SCHEMA; Schema: -; Owner: USERNAMEPLACEHOLDER
--

CREATE SCHEMA fits;


ALTER SCHEMA fits OWNER TO USERNAMEPLACEHOLDER;

--
-- Name: gnn; Type: SCHEMA; Schema: -; Owner: USERNAMEPLACEHOLDER
--

CREATE SCHEMA gnn;


ALTER SCHEMA gnn OWNER TO USERNAMEPLACEHOLDER;

--
-- Name: kdg; Type: SCHEMA; Schema: -; Owner: USERNAMEPLACEHOLDER
--

CREATE SCHEMA kdg;


ALTER SCHEMA kdg OWNER TO USERNAMEPLACEHOLDER;

--
-- Name: rdkit; Type: EXTENSION; Schema: -; Owner: -
--

CREATE EXTENSION IF NOT EXISTS rdkit WITH SCHEMA public;

--
-- Name: EXTENSION rdkit; Type: COMMENT; Schema: -; Owner: 
--

COMMENT ON EXTENSION rdkit IS 'Cheminformatics functionality for PostgreSQL.';


SET default_tablespace = '';

SET default_table_access_method = heap;

--
-- Name: mapping; Type: TABLE; Schema: bookkeeping; Owner: USERNAMEPLACEHOLDER
--

CREATE TABLE bookkeeping.mapping (
    molregno integer NOT NULL,
    idx integer
);


ALTER TABLE bookkeeping.mapping OWNER TO USERNAMEPLACEHOLDER;

--
-- Name: comparisonsmetadata; Type: TABLE; Schema: comparisons; Owner: USERNAMEPLACEHOLDER
--

CREATE TABLE comparisons.comparisonsmetadata (
    comparisonid integer NOT NULL,
    datasourcea text,
    hierarchya boolean,
    datasourceb text,
    hierarchyb boolean,
    environmenta text,
    environmentb text
);


ALTER TABLE comparisons.comparisonsmetadata OWNER TO USERNAMEPLACEHOLDER;

--
-- Name: et4et2; Type: TABLE; Schema: comparisons; Owner: USERNAMEPLACEHOLDER
--

CREATE TABLE comparisons.et4et2 (
    entrynr integer NOT NULL,
    hierarchy boolean,
    multiplicityet2 integer,
    multiplicityet4 integer,
    peakpositionset2 integer[],
    peakpositionset4 integer[],
    solvent text,
    note text
);


ALTER TABLE comparisons.et4et2 OWNER TO USERNAMEPLACEHOLDER;

--
-- Name: fits; Type: TABLE; Schema: comparisons; Owner: USERNAMEPLACEHOLDER
--

CREATE TABLE comparisons.fits (
    entrynr integer NOT NULL,
    env text NOT NULL,
    energythreshold integer NOT NULL,
    coeffs numeric[]
);


ALTER TABLE comparisons.fits OWNER TO USERNAMEPLACEHOLDER;

--
-- Name: metrics; Type: TABLE; Schema: comparisons; Owner: USERNAMEPLACEHOLDER
--

CREATE TABLE comparisons.metrics (
    comparisonid integer NOT NULL,
    entrynr integer NOT NULL,
    threshold integer,
    omega double precision,
    na integer,
    nb integer,
    wasserstein2 double precision,
    wasserstein1 double precision,
    tau double precision,
    omeganorm double precision,
    wasserstein2norm double precision,
    wasserstein1norm double precision
);


ALTER TABLE comparisons.metrics OWNER TO USERNAMEPLACEHOLDER;

--
-- Name: wassersteindistances; Type: TABLE; Schema: comparisons; Owner: USERNAMEPLACEHOLDER
--

CREATE TABLE comparisons.wassersteindistances (
    comparisonid integer NOT NULL,
    entrynr integer NOT NULL,
    wasserstein double precision,
    na integer,
    nb integer
);


ALTER TABLE comparisons.wassersteindistances OWNER TO USERNAMEPLACEHOLDER;

--
-- Name: conformers; Type: TABLE; Schema: csd202403; Owner: USERNAMEPLACEHOLDER
--

CREATE TABLE csd202403.conformers (
    conf_id integer NOT NULL,
    molregno integer,
    conformer_hash text NOT NULL,
    molblock text
);


ALTER TABLE csd202403.conformers OWNER TO USERNAMEPLACEHOLDER;

--
-- Name: conformers_conf_id_seq; Type: SEQUENCE; Schema: csd202403; Owner: USERNAMEPLACEHOLDER
--

CREATE SEQUENCE csd202403.conformers_conf_id_seq
    AS integer
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


ALTER SEQUENCE csd202403.conformers_conf_id_seq OWNER TO USERNAMEPLACEHOLDER;

--
-- Name: conformers_conf_id_seq; Type: SEQUENCE OWNED BY; Schema: csd202403; Owner: USERNAMEPLACEHOLDER
--

ALTER SEQUENCE csd202403.conformers_conf_id_seq OWNED BY csd202403.conformers.conf_id;


--
-- Name: hashes; Type: TABLE; Schema: csd202403; Owner: USERNAMEPLACEHOLDER
--

CREATE TABLE csd202403.hashes (
    molregno integer NOT NULL,
    fullhash text,
    formula text,
    canonical_smiles text,
    no_stereo_smiles text,
    tautomer_hash text,
    no_stereo_tautomer_hash text,
    escape text,
    sgroup_data text,
    rdkitversion text
);


ALTER TABLE csd202403.hashes OWNER TO USERNAMEPLACEHOLDER;

--
-- Name: hashes_molregno_seq; Type: SEQUENCE; Schema: csd202403; Owner: USERNAMEPLACEHOLDER
--

CREATE SEQUENCE csd202403.hashes_molregno_seq
    AS integer
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


ALTER SEQUENCE csd202403.hashes_molregno_seq OWNER TO USERNAMEPLACEHOLDER;

--
-- Name: hashes_molregno_seq; Type: SEQUENCE OWNED BY; Schema: csd202403; Owner: USERNAMEPLACEHOLDER
--

ALTER SEQUENCE csd202403.hashes_molregno_seq OWNED BY csd202403.hashes.molregno;


--
-- Name: molblocks; Type: TABLE; Schema: csd202403; Owner: USERNAMEPLACEHOLDER
--

CREATE TABLE csd202403.molblocks (
    molregno integer NOT NULL,
    molblock text,
    standardization text
);


ALTER TABLE csd202403.molblocks OWNER TO USERNAMEPLACEHOLDER;

--
-- Name: orig_data; Type: TABLE; Schema: csd202403; Owner: USERNAMEPLACEHOLDER
--

CREATE TABLE csd202403.orig_data (
    molregno integer NOT NULL,
    data text,
    datatype text,
    "timestamp" timestamp without time zone DEFAULT now()
);


ALTER TABLE csd202403.orig_data OWNER TO USERNAMEPLACEHOLDER;

--
-- Name: profiles; Type: TABLE; Schema: csd202403; Owner: USERNAMEPLACEHOLDER
--

CREATE TABLE csd202403.profiles (
    confid integer NOT NULL,
    dihedral integer[] NOT NULL,
    torsionvalue numeric
);


ALTER TABLE csd202403.profiles OWNER TO USERNAMEPLACEHOLDER;

--
-- Name: registration_metadata; Type: TABLE; Schema: csd202403; Owner: USERNAMEPLACEHOLDER
--

CREATE TABLE csd202403.registration_metadata (
    key text,
    value text
);


ALTER TABLE csd202403.registration_metadata OWNER TO USERNAMEPLACEHOLDER;

--
-- Name: caches; Type: TABLE; Schema: csd202403modified; Owner: USERNAMEPLACEHOLDER
--

CREATE TABLE csd202403modified.caches (
    molregno integer,
    cache text
);


ALTER TABLE csd202403modified.caches OWNER TO USERNAMEPLACEHOLDER;

--
-- Name: conformers; Type: TABLE; Schema: csd202403modified; Owner: USERNAMEPLACEHOLDER
--

CREATE TABLE csd202403modified.conformers (
    conf_id integer NOT NULL,
    molregno integer,
    conformer_hash text NOT NULL,
    molblock text
);


ALTER TABLE csd202403modified.conformers OWNER TO USERNAMEPLACEHOLDER;

--
-- Name: conformers_conf_id_seq; Type: SEQUENCE; Schema: csd202403modified; Owner: USERNAMEPLACEHOLDER
--

CREATE SEQUENCE csd202403modified.conformers_conf_id_seq
    AS integer
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


ALTER SEQUENCE csd202403modified.conformers_conf_id_seq OWNER TO USERNAMEPLACEHOLDER;

--
-- Name: conformers_conf_id_seq; Type: SEQUENCE OWNED BY; Schema: csd202403modified; Owner: USERNAMEPLACEHOLDER
--

ALTER SEQUENCE csd202403modified.conformers_conf_id_seq OWNED BY csd202403modified.conformers.conf_id;


--
-- Name: csd202403org; Type: TABLE; Schema: csd202403modified; Owner: USERNAMEPLACEHOLDER
--

CREATE TABLE csd202403modified.csd202403org (
    molregnomod integer,
    confidmod integer,
    molregnoorg integer,
    confidorg integer
);


ALTER TABLE csd202403modified.csd202403org OWNER TO USERNAMEPLACEHOLDER;

--
-- Name: hashes; Type: TABLE; Schema: csd202403modified; Owner: USERNAMEPLACEHOLDER
--

CREATE TABLE csd202403modified.hashes (
    molregno integer NOT NULL,
    fullhash text,
    formula text,
    canonical_smiles text,
    no_stereo_smiles text,
    tautomer_hash text,
    no_stereo_tautomer_hash text,
    escape text,
    sgroup_data text,
    rdkitversion text
);


ALTER TABLE csd202403modified.hashes OWNER TO USERNAMEPLACEHOLDER;

--
-- Name: hashes_molregno_seq; Type: SEQUENCE; Schema: csd202403modified; Owner: USERNAMEPLACEHOLDER
--

CREATE SEQUENCE csd202403modified.hashes_molregno_seq
    AS integer
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


ALTER SEQUENCE csd202403modified.hashes_molregno_seq OWNER TO USERNAMEPLACEHOLDER;

--
-- Name: hashes_molregno_seq; Type: SEQUENCE OWNED BY; Schema: csd202403modified; Owner: USERNAMEPLACEHOLDER
--

ALTER SEQUENCE csd202403modified.hashes_molregno_seq OWNED BY csd202403modified.hashes.molregno;


--
-- Name: molblocks; Type: TABLE; Schema: csd202403modified; Owner: USERNAMEPLACEHOLDER
--

CREATE TABLE csd202403modified.molblocks (
    molregno integer NOT NULL,
    molblock text,
    standardization text
);


ALTER TABLE csd202403modified.molblocks OWNER TO USERNAMEPLACEHOLDER;

--
-- Name: orig_data; Type: TABLE; Schema: csd202403modified; Owner: USERNAMEPLACEHOLDER
--

CREATE TABLE csd202403modified.orig_data (
    molregno integer NOT NULL,
    data text,
    datatype text,
    "timestamp" timestamp without time zone DEFAULT now()
);


ALTER TABLE csd202403modified.orig_data OWNER TO USERNAMEPLACEHOLDER;

--
-- Name: profiles; Type: TABLE; Schema: csd202403modified; Owner: USERNAMEPLACEHOLDER
--

CREATE TABLE csd202403modified.profiles (
    molregno integer NOT NULL,
    dihedral integer[] NOT NULL,
    torsionvalueradian numeric,
    torsionvalue numeric
);


ALTER TABLE csd202403modified.profiles OWNER TO USERNAMEPLACEHOLDER;

--
-- Name: registration_metadata; Type: TABLE; Schema: csd202403modified; Owner: USERNAMEPLACEHOLDER
--

CREATE TABLE csd202403modified.registration_metadata (
    key text,
    value text
);


ALTER TABLE csd202403modified.registration_metadata OWNER TO USERNAMEPLACEHOLDER;

--
-- Name: mols; Type: TABLE; Schema: csd202403modifiedrdk; Owner: USERNAMEPLACEHOLDER
--

CREATE TABLE csd202403modifiedrdk.mols (
    molregno integer,
    m public.mol
);


ALTER TABLE csd202403modifiedrdk.mols OWNER TO USERNAMEPLACEHOLDER;

--
-- Name: mols; Type: TABLE; Schema: csd202403rdk; Owner: USERNAMEPLACEHOLDER
--

CREATE TABLE csd202403rdk.mols (
    molregno integer,
    m public.mol
);


ALTER TABLE csd202403rdk.mols OWNER TO USERNAMEPLACEHOLDER;

--
-- Name: et2; Type: TABLE; Schema: extractions; Owner: USERNAMEPLACEHOLDER
--

CREATE TABLE extractions.et2 (
    entrynr integer,
    smarts text,
    molregno integer,
    dihedral integer[],
    hierarchy boolean
);


ALTER TABLE extractions.et2 OWNER TO USERNAMEPLACEHOLDER;

--
-- Name: et4; Type: TABLE; Schema: extractions; Owner: USERNAMEPLACEHOLDER
--

CREATE TABLE extractions.et4 (
    entrynr integer NOT NULL,
    molregno integer NOT NULL,
    dihedral integer[] NOT NULL,
    hierarchy boolean NOT NULL,
    bond integer[]
);


ALTER TABLE extractions.et4 OWNER TO USERNAMEPLACEHOLDER;

--
-- Name: et4csd202403; Type: TABLE; Schema: extractions; Owner: USERNAMEPLACEHOLDER
--

CREATE TABLE extractions.et4csd202403 (
    entrynr integer NOT NULL,
    confid integer NOT NULL,
    dihedral integer[] NOT NULL,
    hierarchy boolean NOT NULL,
    bond integer[]
);


ALTER TABLE extractions.et4csd202403 OWNER TO USERNAMEPLACEHOLDER;

--
-- Name: et4csd202403modified; Type: TABLE; Schema: extractions; Owner: USERNAMEPLACEHOLDER
--

CREATE TABLE extractions.et4csd202403modified (
    entrynr integer NOT NULL,
    molregno integer NOT NULL,
    dihedral integer[] NOT NULL,
    hierarchy boolean NOT NULL,
    bond integer[]
);


ALTER TABLE extractions.et4csd202403modified OWNER TO USERNAMEPLACEHOLDER;

--
-- Name: et4hierarchy; Type: TABLE; Schema: extractions; Owner: USERNAMEPLACEHOLDER
--

CREATE TABLE extractions.et4hierarchy (
    entrynr integer NOT NULL,
    smarts text,
    hierarchyclass text,
    subhierarchyclass text,
    ignored boolean,
    note text,
    alteredsmarts text
);


ALTER TABLE extractions.et4hierarchy OWNER TO USERNAMEPLACEHOLDER;

--
-- Name: fits; Type: TABLE; Schema: fits; Owner: USERNAMEPLACEHOLDER
--

CREATE TABLE fits.fits (
    type text,
    fitparameters numeric[],
    ngaussians integer,
    peakpositions numeric[],
    entrynr integer,
    nbins integer,
    hierarchy boolean,
    etversion integer
);


ALTER TABLE fits.fits OWNER TO USERNAMEPLACEHOLDER;

--
-- Name: chloroform; Type: TABLE; Schema: gnn; Owner: USERNAMEPLACEHOLDER
--

CREATE TABLE gnn.chloroform (
    molregno integer NOT NULL,
    dihedral integer[] NOT NULL,
    torsionvalues numeric[]
);


ALTER TABLE gnn.chloroform OWNER TO USERNAMEPLACEHOLDER;

--
-- Name: chloroformmetadata; Type: TABLE; Schema: gnn; Owner: USERNAMEPLACEHOLDER
--

CREATE TABLE gnn.chloroformmetadata (
    molregno integer NOT NULL,
    coords numeric[],
    energies numeric[],
    energiesnormalized numeric[]
);


ALTER TABLE gnn.chloroformmetadata OWNER TO USERNAMEPLACEHOLDER;

--
-- Name: csdhexane; Type: TABLE; Schema: gnn; Owner: USERNAMEPLACEHOLDER
--

CREATE TABLE gnn.csdhexane (
    molregno integer NOT NULL,
    dihedral integer[] NOT NULL,
    torsionvalues numeric[]
);


ALTER TABLE gnn.csdhexane OWNER TO USERNAMEPLACEHOLDER;

--
-- Name: csdhexanemetadata; Type: TABLE; Schema: gnn; Owner: USERNAMEPLACEHOLDER
--

CREATE TABLE gnn.csdhexanemetadata (
    molregno integer NOT NULL,
    coords numeric[],
    energies numeric[],
    energiesnormalized numeric[],
    error boolean
);


ALTER TABLE gnn.csdhexanemetadata OWNER TO USERNAMEPLACEHOLDER;

--
-- Name: csdtip3p; Type: TABLE; Schema: gnn; Owner: USERNAMEPLACEHOLDER
--

CREATE TABLE gnn.csdtip3p (
    molregno integer NOT NULL,
    dihedral integer[] NOT NULL,
    torsionvalues numeric[]
);


ALTER TABLE gnn.csdtip3p OWNER TO USERNAMEPLACEHOLDER;

--
-- Name: csdtip3pmetadata; Type: TABLE; Schema: gnn; Owner: USERNAMEPLACEHOLDER
--

CREATE TABLE gnn.csdtip3pmetadata (
    molregno integer NOT NULL,
    coords numeric[],
    energies numeric[],
    energiesnormalized numeric[],
    error boolean
);


ALTER TABLE gnn.csdtip3pmetadata OWNER TO USERNAMEPLACEHOLDER;

--
-- Name: csdvac; Type: TABLE; Schema: gnn; Owner: USERNAMEPLACEHOLDER
--

CREATE TABLE gnn.csdvac (
    molregno integer NOT NULL,
    dihedral integer[] NOT NULL,
    torsionvalues numeric[]
);


ALTER TABLE gnn.csdvac OWNER TO USERNAMEPLACEHOLDER;

--
-- Name: csdvacmetadata; Type: TABLE; Schema: gnn; Owner: USERNAMEPLACEHOLDER
--

CREATE TABLE gnn.csdvacmetadata (
    molregno integer NOT NULL,
    coords numeric[],
    energies numeric[],
    energiesnormalized numeric[],
    error boolean
);


ALTER TABLE gnn.csdvacmetadata OWNER TO USERNAMEPLACEHOLDER;

--
-- Name: dmso; Type: TABLE; Schema: gnn; Owner: USERNAMEPLACEHOLDER
--

CREATE TABLE gnn.dmso (
    molregno integer NOT NULL,
    dihedral integer[] NOT NULL,
    torsionvalues numeric[]
);


ALTER TABLE gnn.dmso OWNER TO USERNAMEPLACEHOLDER;

--
-- Name: dmsometadata; Type: TABLE; Schema: gnn; Owner: USERNAMEPLACEHOLDER
--

CREATE TABLE gnn.dmsometadata (
    molregno integer NOT NULL,
    coords numeric[],
    energies numeric[],
    energiesnormalized numeric[]
);


ALTER TABLE gnn.dmsometadata OWNER TO USERNAMEPLACEHOLDER;

--
-- Name: hexane; Type: TABLE; Schema: gnn; Owner: USERNAMEPLACEHOLDER
--

CREATE TABLE gnn.hexane (
    molregno integer NOT NULL,
    dihedral integer[] NOT NULL,
    torsionvalues numeric[]
);


ALTER TABLE gnn.hexane OWNER TO USERNAMEPLACEHOLDER;

--
-- Name: hexanemetadata; Type: TABLE; Schema: gnn; Owner: USERNAMEPLACEHOLDER
--

CREATE TABLE gnn.hexanemetadata (
    molregno integer NOT NULL,
    coords numeric[],
    energies numeric[],
    energiesnormalized numeric[]
);


ALTER TABLE gnn.hexanemetadata OWNER TO USERNAMEPLACEHOLDER;

--
-- Name: methanol; Type: TABLE; Schema: gnn; Owner: USERNAMEPLACEHOLDER
--

CREATE TABLE gnn.methanol (
    molregno integer NOT NULL,
    dihedral integer[] NOT NULL,
    torsionvalues numeric[]
);


ALTER TABLE gnn.methanol OWNER TO USERNAMEPLACEHOLDER;

--
-- Name: methanolmetadata; Type: TABLE; Schema: gnn; Owner: USERNAMEPLACEHOLDER
--

CREATE TABLE gnn.methanolmetadata (
    molregno integer NOT NULL,
    coords numeric[],
    energies numeric[],
    energiesnormalized numeric[]
);


ALTER TABLE gnn.methanolmetadata OWNER TO USERNAMEPLACEHOLDER;

--
-- Name: tip3p; Type: TABLE; Schema: gnn; Owner: USERNAMEPLACEHOLDER
--

CREATE TABLE gnn.tip3p (
    molregno integer NOT NULL,
    dihedral integer[] NOT NULL,
    torsionvalues numeric[]
);


ALTER TABLE gnn.tip3p OWNER TO USERNAMEPLACEHOLDER;

--
-- Name: tip3pmetadata; Type: TABLE; Schema: gnn; Owner: USERNAMEPLACEHOLDER
--

CREATE TABLE gnn.tip3pmetadata (
    molregno integer NOT NULL,
    coords numeric[],
    energies numeric[],
    energiesnormalized numeric[]
);


ALTER TABLE gnn.tip3pmetadata OWNER TO USERNAMEPLACEHOLDER;

--
-- Name: vac; Type: TABLE; Schema: gnn; Owner: USERNAMEPLACEHOLDER
--

CREATE TABLE gnn.vac (
    molregno integer NOT NULL,
    dihedral integer[] NOT NULL,
    torsionvalues numeric[]
);


ALTER TABLE gnn.vac OWNER TO USERNAMEPLACEHOLDER;

--
-- Name: vacmetadata; Type: TABLE; Schema: gnn; Owner: USERNAMEPLACEHOLDER
--

CREATE TABLE gnn.vacmetadata (
    molregno integer NOT NULL,
    coords numeric[],
    energies numeric[],
    energiesnormalized numeric[],
    error boolean
);


ALTER TABLE gnn.vacmetadata OWNER TO USERNAMEPLACEHOLDER;

--
-- Name: metadata; Type: TABLE; Schema: kdg; Owner: USERNAMEPLACEHOLDER
--

CREATE TABLE kdg.metadata (
    molregno integer NOT NULL,
    coords numeric[]
);


ALTER TABLE kdg.metadata OWNER TO USERNAMEPLACEHOLDER;

--
-- Name: profiles; Type: TABLE; Schema: kdg; Owner: USERNAMEPLACEHOLDER
--

CREATE TABLE kdg.profiles (
    molregno integer NOT NULL,
    dihedral integer[] NOT NULL,
    torsionvalues numeric[]
);


ALTER TABLE kdg.profiles OWNER TO USERNAMEPLACEHOLDER;

--
-- Name: hashes; Type: TABLE; Schema: public; Owner: USERNAMEPLACEHOLDER
--

CREATE TABLE public.hashes (
    molregno integer NOT NULL,
    fullhash text,
    formula text,
    canonical_smiles text,
    no_stereo_smiles text,
    tautomer_hash text,
    no_stereo_tautomer_hash text,
    escape text,
    sgroup_data text,
    rdkitversion text
);


ALTER TABLE public.hashes OWNER TO USERNAMEPLACEHOLDER;

--
-- Name: hashes_molregno_seq; Type: SEQUENCE; Schema: public; Owner: USERNAMEPLACEHOLDER
--

CREATE SEQUENCE public.hashes_molregno_seq
    AS integer
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


ALTER SEQUENCE public.hashes_molregno_seq OWNER TO USERNAMEPLACEHOLDER;

--
-- Name: hashes_molregno_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: USERNAMEPLACEHOLDER
--

ALTER SEQUENCE public.hashes_molregno_seq OWNED BY public.hashes.molregno;


--
-- Name: molblocks; Type: TABLE; Schema: public; Owner: USERNAMEPLACEHOLDER
--

CREATE TABLE public.molblocks (
    molregno integer NOT NULL,
    molblock text,
    standardization text
);


ALTER TABLE public.molblocks OWNER TO USERNAMEPLACEHOLDER;

--
-- Name: orig_data; Type: TABLE; Schema: public; Owner: USERNAMEPLACEHOLDER
--

CREATE TABLE public.orig_data (
    molregno integer NOT NULL,
    data text,
    datatype text,
    "timestamp" timestamp without time zone DEFAULT now()
);


ALTER TABLE public.orig_data OWNER TO USERNAMEPLACEHOLDER;

--
-- Name: registration_metadata; Type: TABLE; Schema: public; Owner: USERNAMEPLACEHOLDER
--

CREATE TABLE public.registration_metadata (
    key text,
    value text
);


ALTER TABLE public.registration_metadata OWNER TO USERNAMEPLACEHOLDER;

--
-- Name: conformers conf_id; Type: DEFAULT; Schema: csd202403; Owner: USERNAMEPLACEHOLDER
--

ALTER TABLE ONLY csd202403.conformers ALTER COLUMN conf_id SET DEFAULT nextval('csd202403.conformers_conf_id_seq'::regclass);


--
-- Name: hashes molregno; Type: DEFAULT; Schema: csd202403; Owner: USERNAMEPLACEHOLDER
--

ALTER TABLE ONLY csd202403.hashes ALTER COLUMN molregno SET DEFAULT nextval('csd202403.hashes_molregno_seq'::regclass);


--
-- Name: conformers conf_id; Type: DEFAULT; Schema: csd202403modified; Owner: USERNAMEPLACEHOLDER
--

ALTER TABLE ONLY csd202403modified.conformers ALTER COLUMN conf_id SET DEFAULT nextval('csd202403modified.conformers_conf_id_seq'::regclass);


--
-- Name: hashes molregno; Type: DEFAULT; Schema: csd202403modified; Owner: USERNAMEPLACEHOLDER
--

ALTER TABLE ONLY csd202403modified.hashes ALTER COLUMN molregno SET DEFAULT nextval('csd202403modified.hashes_molregno_seq'::regclass);


--
-- Name: hashes molregno; Type: DEFAULT; Schema: public; Owner: USERNAMEPLACEHOLDER
--

ALTER TABLE ONLY public.hashes ALTER COLUMN molregno SET DEFAULT nextval('public.hashes_molregno_seq'::regclass);


--
-- Name: mapping mapping_pkey; Type: CONSTRAINT; Schema: bookkeeping; Owner: USERNAMEPLACEHOLDER
--

ALTER TABLE ONLY bookkeeping.mapping
    ADD CONSTRAINT mapping_pkey PRIMARY KEY (molregno);


--
-- Name: comparisonsmetadata comparisonsmetadata_pkey; Type: CONSTRAINT; Schema: comparisons; Owner: USERNAMEPLACEHOLDER
--

ALTER TABLE ONLY comparisons.comparisonsmetadata
    ADD CONSTRAINT comparisonsmetadata_pkey PRIMARY KEY (comparisonid);


--
-- Name: et4et2 et4et2_pkey; Type: CONSTRAINT; Schema: comparisons; Owner: USERNAMEPLACEHOLDER
--

ALTER TABLE ONLY comparisons.et4et2
    ADD CONSTRAINT et4et2_pkey PRIMARY KEY (entrynr);


--
-- Name: fits fits_pk; Type: CONSTRAINT; Schema: comparisons; Owner: USERNAMEPLACEHOLDER
--

ALTER TABLE ONLY comparisons.fits
    ADD CONSTRAINT fits_pk PRIMARY KEY (entrynr, env, energythreshold);


--
-- Name: metrics metrics_pkey; Type: CONSTRAINT; Schema: comparisons; Owner: USERNAMEPLACEHOLDER
--

ALTER TABLE ONLY comparisons.metrics
    ADD CONSTRAINT metrics_pkey PRIMARY KEY (comparisonid, entrynr);


--
-- Name: wassersteindistances wassersteindistances_pkey; Type: CONSTRAINT; Schema: comparisons; Owner: USERNAMEPLACEHOLDER
--

ALTER TABLE ONLY comparisons.wassersteindistances
    ADD CONSTRAINT wassersteindistances_pkey PRIMARY KEY (comparisonid, entrynr);


--
-- Name: conformers conformers_conformer_hash_key; Type: CONSTRAINT; Schema: csd202403; Owner: USERNAMEPLACEHOLDER
--

ALTER TABLE ONLY csd202403.conformers
    ADD CONSTRAINT conformers_conformer_hash_key UNIQUE (conformer_hash);


--
-- Name: conformers conformers_pkey; Type: CONSTRAINT; Schema: csd202403; Owner: USERNAMEPLACEHOLDER
--

ALTER TABLE ONLY csd202403.conformers
    ADD CONSTRAINT conformers_pkey PRIMARY KEY (conf_id);


--
-- Name: hashes hashes_fullhash_key; Type: CONSTRAINT; Schema: csd202403; Owner: USERNAMEPLACEHOLDER
--

ALTER TABLE ONLY csd202403.hashes
    ADD CONSTRAINT hashes_fullhash_key UNIQUE (fullhash);


--
-- Name: hashes hashes_pkey; Type: CONSTRAINT; Schema: csd202403; Owner: USERNAMEPLACEHOLDER
--

ALTER TABLE ONLY csd202403.hashes
    ADD CONSTRAINT hashes_pkey PRIMARY KEY (molregno);


--
-- Name: molblocks molblocks_molregno_key; Type: CONSTRAINT; Schema: csd202403; Owner: USERNAMEPLACEHOLDER
--

ALTER TABLE ONLY csd202403.molblocks
    ADD CONSTRAINT molblocks_molregno_key UNIQUE (molregno);


--
-- Name: orig_data orig_data_molregno_key; Type: CONSTRAINT; Schema: csd202403; Owner: USERNAMEPLACEHOLDER
--

ALTER TABLE ONLY csd202403.orig_data
    ADD CONSTRAINT orig_data_molregno_key UNIQUE (molregno);


--
-- Name: profiles profiles_pkey; Type: CONSTRAINT; Schema: csd202403; Owner: USERNAMEPLACEHOLDER
--

ALTER TABLE ONLY csd202403.profiles
    ADD CONSTRAINT profiles_pkey PRIMARY KEY (confid, dihedral);


--
-- Name: conformers conformers_conformer_hash_key; Type: CONSTRAINT; Schema: csd202403modified; Owner: USERNAMEPLACEHOLDER
--

ALTER TABLE ONLY csd202403modified.conformers
    ADD CONSTRAINT conformers_conformer_hash_key UNIQUE (conformer_hash);


--
-- Name: conformers conformers_pkey; Type: CONSTRAINT; Schema: csd202403modified; Owner: USERNAMEPLACEHOLDER
--

ALTER TABLE ONLY csd202403modified.conformers
    ADD CONSTRAINT conformers_pkey PRIMARY KEY (conf_id);


--
-- Name: hashes hashes_fullhash_key; Type: CONSTRAINT; Schema: csd202403modified; Owner: USERNAMEPLACEHOLDER
--

ALTER TABLE ONLY csd202403modified.hashes
    ADD CONSTRAINT hashes_fullhash_key UNIQUE (fullhash);


--
-- Name: hashes hashes_pkey; Type: CONSTRAINT; Schema: csd202403modified; Owner: USERNAMEPLACEHOLDER
--

ALTER TABLE ONLY csd202403modified.hashes
    ADD CONSTRAINT hashes_pkey PRIMARY KEY (molregno);


--
-- Name: molblocks molblocks_molregno_key; Type: CONSTRAINT; Schema: csd202403modified; Owner: USERNAMEPLACEHOLDER
--

ALTER TABLE ONLY csd202403modified.molblocks
    ADD CONSTRAINT molblocks_molregno_key UNIQUE (molregno);


--
-- Name: orig_data orig_data_molregno_key; Type: CONSTRAINT; Schema: csd202403modified; Owner: USERNAMEPLACEHOLDER
--

ALTER TABLE ONLY csd202403modified.orig_data
    ADD CONSTRAINT orig_data_molregno_key UNIQUE (molregno);


--
-- Name: profiles profiles_pkey; Type: CONSTRAINT; Schema: csd202403modified; Owner: USERNAMEPLACEHOLDER
--

ALTER TABLE ONLY csd202403modified.profiles
    ADD CONSTRAINT profiles_pkey PRIMARY KEY (molregno, dihedral);


--
-- Name: et4hierarchy et4_pkey; Type: CONSTRAINT; Schema: extractions; Owner: USERNAMEPLACEHOLDER
--

ALTER TABLE ONLY extractions.et4hierarchy
    ADD CONSTRAINT et4_pkey PRIMARY KEY (entrynr);


--
-- Name: et4 et4_pkey1; Type: CONSTRAINT; Schema: extractions; Owner: USERNAMEPLACEHOLDER
--

ALTER TABLE ONLY extractions.et4
    ADD CONSTRAINT et4_pkey1 PRIMARY KEY (entrynr, molregno, dihedral, hierarchy);


--
-- Name: et4csd202403 et4csd202403_pkey; Type: CONSTRAINT; Schema: extractions; Owner: USERNAMEPLACEHOLDER
--

ALTER TABLE ONLY extractions.et4csd202403
    ADD CONSTRAINT et4csd202403_pkey PRIMARY KEY (entrynr, confid, dihedral, hierarchy);


--
-- Name: et4csd202403modified et4csd202403modified_pkey; Type: CONSTRAINT; Schema: extractions; Owner: USERNAMEPLACEHOLDER
--

ALTER TABLE ONLY extractions.et4csd202403modified
    ADD CONSTRAINT et4csd202403modified_pkey PRIMARY KEY (entrynr, molregno, dihedral, hierarchy);


--
-- Name: chloroform chloroform_pkey; Type: CONSTRAINT; Schema: gnn; Owner: USERNAMEPLACEHOLDER
--

ALTER TABLE ONLY gnn.chloroform
    ADD CONSTRAINT chloroform_pkey PRIMARY KEY (molregno, dihedral);


--
-- Name: chloroformmetadata chloroformmetadata_pkey; Type: CONSTRAINT; Schema: gnn; Owner: USERNAMEPLACEHOLDER
--

ALTER TABLE ONLY gnn.chloroformmetadata
    ADD CONSTRAINT chloroformmetadata_pkey PRIMARY KEY (molregno);


--
-- Name: csdhexane csdhexane_pkey; Type: CONSTRAINT; Schema: gnn; Owner: USERNAMEPLACEHOLDER
--

ALTER TABLE ONLY gnn.csdhexane
    ADD CONSTRAINT csdhexane_pkey PRIMARY KEY (molregno, dihedral);


--
-- Name: csdhexanemetadata csdhexanemetadata_pkey; Type: CONSTRAINT; Schema: gnn; Owner: USERNAMEPLACEHOLDER
--

ALTER TABLE ONLY gnn.csdhexanemetadata
    ADD CONSTRAINT csdhexanemetadata_pkey PRIMARY KEY (molregno);


--
-- Name: csdtip3p csdtip3p_pkey; Type: CONSTRAINT; Schema: gnn; Owner: USERNAMEPLACEHOLDER
--

ALTER TABLE ONLY gnn.csdtip3p
    ADD CONSTRAINT csdtip3p_pkey PRIMARY KEY (molregno, dihedral);


--
-- Name: csdtip3pmetadata csdtip3pmetadata_pkey; Type: CONSTRAINT; Schema: gnn; Owner: USERNAMEPLACEHOLDER
--

ALTER TABLE ONLY gnn.csdtip3pmetadata
    ADD CONSTRAINT csdtip3pmetadata_pkey PRIMARY KEY (molregno);


--
-- Name: csdvac csdvac_pkey; Type: CONSTRAINT; Schema: gnn; Owner: USERNAMEPLACEHOLDER
--

ALTER TABLE ONLY gnn.csdvac
    ADD CONSTRAINT csdvac_pkey PRIMARY KEY (molregno, dihedral);


--
-- Name: csdvacmetadata csdvacmetadata_pkey; Type: CONSTRAINT; Schema: gnn; Owner: USERNAMEPLACEHOLDER
--

ALTER TABLE ONLY gnn.csdvacmetadata
    ADD CONSTRAINT csdvacmetadata_pkey PRIMARY KEY (molregno);


--
-- Name: dmso dmso_pkey; Type: CONSTRAINT; Schema: gnn; Owner: USERNAMEPLACEHOLDER
--

ALTER TABLE ONLY gnn.dmso
    ADD CONSTRAINT dmso_pkey PRIMARY KEY (molregno, dihedral);


--
-- Name: dmsometadata dmsometadata_pkey; Type: CONSTRAINT; Schema: gnn; Owner: USERNAMEPLACEHOLDER
--

ALTER TABLE ONLY gnn.dmsometadata
    ADD CONSTRAINT dmsometadata_pkey PRIMARY KEY (molregno);


--
-- Name: hexane hexane_pkey; Type: CONSTRAINT; Schema: gnn; Owner: USERNAMEPLACEHOLDER
--

ALTER TABLE ONLY gnn.hexane
    ADD CONSTRAINT hexane_pkey PRIMARY KEY (molregno, dihedral);


--
-- Name: hexanemetadata hexanemetadata_pkey; Type: CONSTRAINT; Schema: gnn; Owner: USERNAMEPLACEHOLDER
--

ALTER TABLE ONLY gnn.hexanemetadata
    ADD CONSTRAINT hexanemetadata_pkey PRIMARY KEY (molregno);


--
-- Name: methanol methanol_pkey; Type: CONSTRAINT; Schema: gnn; Owner: USERNAMEPLACEHOLDER
--

ALTER TABLE ONLY gnn.methanol
    ADD CONSTRAINT methanol_pkey PRIMARY KEY (molregno, dihedral);


--
-- Name: methanolmetadata methanolmetadata_pkey; Type: CONSTRAINT; Schema: gnn; Owner: USERNAMEPLACEHOLDER
--

ALTER TABLE ONLY gnn.methanolmetadata
    ADD CONSTRAINT methanolmetadata_pkey PRIMARY KEY (molregno);


--
-- Name: tip3p tip3p_pkey; Type: CONSTRAINT; Schema: gnn; Owner: USERNAMEPLACEHOLDER
--

ALTER TABLE ONLY gnn.tip3p
    ADD CONSTRAINT tip3p_pkey PRIMARY KEY (molregno, dihedral);


--
-- Name: tip3pmetadata tip3pmetadata_pkey; Type: CONSTRAINT; Schema: gnn; Owner: USERNAMEPLACEHOLDER
--

ALTER TABLE ONLY gnn.tip3pmetadata
    ADD CONSTRAINT tip3pmetadata_pkey PRIMARY KEY (molregno);


--
-- Name: vac vac_pkey; Type: CONSTRAINT; Schema: gnn; Owner: USERNAMEPLACEHOLDER
--

ALTER TABLE ONLY gnn.vac
    ADD CONSTRAINT vac_pkey PRIMARY KEY (molregno, dihedral);


--
-- Name: vacmetadata vacmetadata_pkey; Type: CONSTRAINT; Schema: gnn; Owner: USERNAMEPLACEHOLDER
--

ALTER TABLE ONLY gnn.vacmetadata
    ADD CONSTRAINT vacmetadata_pkey PRIMARY KEY (molregno);


--
-- Name: metadata metadata_pkey; Type: CONSTRAINT; Schema: kdg; Owner: USERNAMEPLACEHOLDER
--

ALTER TABLE ONLY kdg.metadata
    ADD CONSTRAINT metadata_pkey PRIMARY KEY (molregno);


--
-- Name: profiles profiles_pkey; Type: CONSTRAINT; Schema: kdg; Owner: USERNAMEPLACEHOLDER
--

ALTER TABLE ONLY kdg.profiles
    ADD CONSTRAINT profiles_pkey PRIMARY KEY (molregno, dihedral);


--
-- Name: hashes hashes_fullhash_key; Type: CONSTRAINT; Schema: public; Owner: USERNAMEPLACEHOLDER
--

ALTER TABLE ONLY public.hashes
    ADD CONSTRAINT hashes_fullhash_key UNIQUE (fullhash);


--
-- Name: hashes hashes_pkey; Type: CONSTRAINT; Schema: public; Owner: USERNAMEPLACEHOLDER
--

ALTER TABLE ONLY public.hashes
    ADD CONSTRAINT hashes_pkey PRIMARY KEY (molregno);


--
-- Name: molblocks molblocks_molregno_key; Type: CONSTRAINT; Schema: public; Owner: USERNAMEPLACEHOLDER
--

ALTER TABLE ONLY public.molblocks
    ADD CONSTRAINT molblocks_molregno_key UNIQUE (molregno);


--
-- Name: orig_data orig_data_molregno_key; Type: CONSTRAINT; Schema: public; Owner: USERNAMEPLACEHOLDER
--

ALTER TABLE ONLY public.orig_data
    ADD CONSTRAINT orig_data_molregno_key UNIQUE (molregno);


--
-- Name: csd202403_conformers_hash_idx; Type: INDEX; Schema: csd202403; Owner: USERNAMEPLACEHOLDER
--

CREATE INDEX csd202403_conformers_hash_idx ON csd202403.conformers USING hash (conformer_hash);


--
-- Name: csd202403_hashes_fullhash_idx; Type: INDEX; Schema: csd202403; Owner: USERNAMEPLACEHOLDER
--

CREATE INDEX csd202403_hashes_fullhash_idx ON csd202403.hashes USING hash (fullhash);


--
-- Name: csd202403modified_conformers_hash_idx; Type: INDEX; Schema: csd202403modified; Owner: USERNAMEPLACEHOLDER
--

CREATE INDEX csd202403modified_conformers_hash_idx ON csd202403modified.conformers USING hash (conformer_hash);


--
-- Name: csd202403modified_hashes_fullhash_idx; Type: INDEX; Schema: csd202403modified; Owner: USERNAMEPLACEHOLDER
--

CREATE INDEX csd202403modified_hashes_fullhash_idx ON csd202403modified.hashes USING hash (fullhash);


--
-- Name: molidx; Type: INDEX; Schema: csd202403modifiedrdk; Owner: USERNAMEPLACEHOLDER
--

CREATE INDEX molidx ON csd202403modifiedrdk.mols USING gist (m);


--
-- Name: molidx; Type: INDEX; Schema: csd202403rdk; Owner: USERNAMEPLACEHOLDER
--

CREATE INDEX molidx ON csd202403rdk.mols USING gist (m);


--
-- Name: test; Type: INDEX; Schema: extractions; Owner: USERNAMEPLACEHOLDER
--

CREATE INDEX test ON extractions.et4 USING btree (molregno);


--
-- Name: test2; Type: INDEX; Schema: extractions; Owner: USERNAMEPLACEHOLDER
--

CREATE INDEX test2 ON extractions.et4 USING btree (bond);


--
-- Name: hashes_fullhash_idx; Type: INDEX; Schema: public; Owner: USERNAMEPLACEHOLDER
--

CREATE INDEX hashes_fullhash_idx ON public.hashes USING hash (fullhash);


--
-- Name: conformers conformers_molregno_fkey; Type: FK CONSTRAINT; Schema: csd202403; Owner: USERNAMEPLACEHOLDER
--

ALTER TABLE ONLY csd202403.conformers
    ADD CONSTRAINT conformers_molregno_fkey FOREIGN KEY (molregno) REFERENCES csd202403.hashes(molregno);


--
-- Name: molblocks molblocks_molregno_fkey; Type: FK CONSTRAINT; Schema: csd202403; Owner: USERNAMEPLACEHOLDER
--

ALTER TABLE ONLY csd202403.molblocks
    ADD CONSTRAINT molblocks_molregno_fkey FOREIGN KEY (molregno) REFERENCES csd202403.hashes(molregno);


--
-- Name: molblocks molblocks_molregno_fkey1; Type: FK CONSTRAINT; Schema: csd202403; Owner: USERNAMEPLACEHOLDER
--

ALTER TABLE ONLY csd202403.molblocks
    ADD CONSTRAINT molblocks_molregno_fkey1 FOREIGN KEY (molregno) REFERENCES csd202403.hashes(molregno);


--
-- Name: orig_data orig_data_molregno_fkey; Type: FK CONSTRAINT; Schema: csd202403; Owner: USERNAMEPLACEHOLDER
--

ALTER TABLE ONLY csd202403.orig_data
    ADD CONSTRAINT orig_data_molregno_fkey FOREIGN KEY (molregno) REFERENCES csd202403.hashes(molregno);


--
-- Name: conformers conformers_molregno_fkey; Type: FK CONSTRAINT; Schema: csd202403modified; Owner: USERNAMEPLACEHOLDER
--

ALTER TABLE ONLY csd202403modified.conformers
    ADD CONSTRAINT conformers_molregno_fkey FOREIGN KEY (molregno) REFERENCES csd202403modified.hashes(molregno);


--
-- Name: molblocks molblocks_molregno_fkey; Type: FK CONSTRAINT; Schema: csd202403modified; Owner: USERNAMEPLACEHOLDER
--

ALTER TABLE ONLY csd202403modified.molblocks
    ADD CONSTRAINT molblocks_molregno_fkey FOREIGN KEY (molregno) REFERENCES csd202403modified.hashes(molregno);


--
-- Name: molblocks molblocks_molregno_fkey1; Type: FK CONSTRAINT; Schema: csd202403modified; Owner: USERNAMEPLACEHOLDER
--

ALTER TABLE ONLY csd202403modified.molblocks
    ADD CONSTRAINT molblocks_molregno_fkey1 FOREIGN KEY (molregno) REFERENCES csd202403modified.hashes(molregno);


--
-- Name: orig_data orig_data_molregno_fkey; Type: FK CONSTRAINT; Schema: csd202403modified; Owner: USERNAMEPLACEHOLDER
--

ALTER TABLE ONLY csd202403modified.orig_data
    ADD CONSTRAINT orig_data_molregno_fkey FOREIGN KEY (molregno) REFERENCES csd202403modified.hashes(molregno);


--
-- Name: molblocks molblocks_molregno_fkey; Type: FK CONSTRAINT; Schema: public; Owner: USERNAMEPLACEHOLDER
--

ALTER TABLE ONLY public.molblocks
    ADD CONSTRAINT molblocks_molregno_fkey FOREIGN KEY (molregno) REFERENCES public.hashes(molregno);


--
-- Name: molblocks molblocks_molregno_fkey1; Type: FK CONSTRAINT; Schema: public; Owner: USERNAMEPLACEHOLDER
--

ALTER TABLE ONLY public.molblocks
    ADD CONSTRAINT molblocks_molregno_fkey1 FOREIGN KEY (molregno) REFERENCES public.hashes(molregno);


--
-- Name: orig_data orig_data_molregno_fkey; Type: FK CONSTRAINT; Schema: public; Owner: USERNAMEPLACEHOLDER
--

ALTER TABLE ONLY public.orig_data
    ADD CONSTRAINT orig_data_molregno_fkey FOREIGN KEY (molregno) REFERENCES public.hashes(molregno);

--
-- PostgreSQL database dump complete
--
