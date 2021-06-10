
CREATE TABLE ganesha_measurements(
id INTEGER PRIMARY KEY,
measure_date TEXT,
sample TEXT,
comment TEXT);

CREATE TABLE protein_ps_samples(
id INTEGER PRIMARY KEY,
creation_date TEXT,
protein_id TEXT,
polysaccharide_id TEXT,
buffer_id INTEGER,
comment TEXT);

CREATE TABLE proteins(
CAS_number TEXT PRIMARY KEY,
name TEXT,
short_name TEXT,
EC_number TEXT,
molecular_weight REAL,
pI REAL,
reference_size REAL,
comment TEXT);

CREATE TABLE polysaccharides(
CAS_number TEXT PRIMARY KEY,
name TEXT,
short_name TEXT,
EC_number TEXT,
molecular_weight REAL,
reference_size REAL,
comment TEXT);

CREATE TABLE buffers(
id INTEGER PRIMARY KEY,
name TEXT,
pH REAL,
salt_name TEXT,
salt_concentration REAL,
visc REAL,
n REAL,
ganesha_id INTEGER,
comment TEXT);



