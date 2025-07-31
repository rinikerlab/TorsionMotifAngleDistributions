--
-- PostgreSQL database dump
--

-- Dumped from database version 16.9 (Ubuntu 16.9-1.pgdg22.04+1)
-- Dumped by pg_dump version 17.5 (Ubuntu 17.5-1.pgdg22.04+1)

SET statement_timeout = 0;
SET lock_timeout = 0;
SET idle_in_transaction_session_timeout = 0;
SET transaction_timeout = 0;
SET client_encoding = 'UTF8';
SET standard_conforming_strings = on;
SELECT pg_catalog.set_config('search_path', '', false);
SET check_function_bodies = false;
SET xmloption = content;
SET client_min_messages = warning;
SET row_security = off;

SET default_tablespace = '';

SET default_table_access_method = heap;

--
-- Data for Name: comparisonsmetadata; Type: TABLE DATA; Schema: comparisons; Owner:
--

COPY comparisons.comparisonsmetadata (comparisonid, datasourcea, hierarchya, datasourceb, hierarchyb, environmenta, environmentb) FROM stdin;
1	csd202403modified	t	csd202403modified	t	crystal	vac
2	csd202403modified	t	csd202403modified	t	crystal	tip3p
3	csd202403modified	t	csd202403modified	t	vac	tip3p
4	csd202403modified	t	csd202403modified	t	crystal	hexane
5	csd202403modified	t	csd202403modified	t	vac	hexane
6	csd202403modified	t	csd202403modified	t	tip3p	hexane
\.

--
-- PostgreSQL database dump complete
--

