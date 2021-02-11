"""
SampleModule example
====================

nananana NaNaNaNa EEEEYOOOO Gooodbye.
"""
#!/usr/bin/env python
# coding: utf-8

# # Borehole Database  
# 
# 
# It is no surprise that the core of a data-driven project is data. Data which is organized in a specific manner, so that users can easily access, analyse, and manipulate data. Many different schemes exist, in which data can be organized, with a most common one being spreadsheet-like tables. Thus, spreadsheet software like Microsoft Excel or Libre-Office Calc are among popular solutions, when it comes to working with data.
# 
# With growing amount of data, however, these software solutions may soon meet their limits, as they can get overly complicated. One example would be many `.xls` files, which are connected among each other using hyperlinks. This is obviously an error-prone solution, not really practical. Thus, greater amounts of data with a more complex structure, are usually maintained in a [database](https://en.wikipedia.org/wiki/Database), following a certain [data model](https://en.wikipedia.org/wiki/Data_model). Here, we use [SQLITE](https://www.sqlite.org/index.html) \cite{hipp_sqlite_2019} as underlying database solution, a SQL database engine. 
# 
# In the following, we will briefly describe the database structure, its content and provide some short examples how to access the database and work with the stored data.

# ## Data model
# 
# Within the database, we follow, as we use SQL, a [relational model](https://en.wikipedia.org/wiki/Relational_model) to organize stored data. This data comprises mainly borehole temperature measurements in the study area. The data was originally compiled by Schärli and Kohl \cite{scharli_archivierung_2002} in a set of excel tables. This *original* data, i.e. in its excel form, is available as supplementary material to the NAGRA Working report [NAB 12-61](https://www.nagra.ch/de/cat/publikationen/arbeitsberichte-nabs/nabs-2012/downloadcenter.htm). This report comprises temperature measurements for boreholes all over Switzerland. Additionally, a stratigraphical description is available for some boreholes. Figure \ref{fig:borehole_map} shows boreholes in Switzerland, which are deeper than 500 m. 
# 
# 
# Many of the temperature data from these deep boreholes is compiled in Schärli and Kohl \cite{scharli_archivierung_2002}, in addition to temperature data from *shallow* boreholes, i.e. shallower than 500 m.
# In this work, we use a subset of this data which is (**a**) inside our area of interest, and (**b**) publicly available data. For instance, figure \ref{fig:database_map} shows a subset of deep boreholes (triangles) in the study area, colored by the data restriction. While blue represents open data, boreholes colored in red contain confidential data. Within the database, this information is stored, so confidential data can easily be erased from the database, in case it is made public.
# 
# 
# Currently, the database contains three related tables:  
# * general borehole information (coordinates, name, original source, ...)  
# * temperature depth information for all boreholes  
# * available petrophysical information 
# 
# Not much petrophysical data is available from the boreholes. Temperature depth information, however, comprises more than 39000 data points. In the following, we present methods and procedures to access these data and work with it from within this notebook. For this, we use a mixture of SQL queries and methods of the data analysis library [`pandas`](https://pandas.pydata.org/). 

# ## Acessing data and visualizing 
# Querying a database is maybe the most often performed task, when it comes to databases. When you type something in a seach bar, for example, you query a database for the words you are looking for. The same, though in a more rudimentary form, can be done with the compiled "borehole temperature" database. 
# 
# The following code cells in this notebook show how:
# * to connect to the database  
# * introduces a very small library `db_access`
# * get information about available tables in the database
# * formulate queries to get desired data, e.g. temperature depth pairs for a specific borehole
# * store query results in a pandas dataframe and visualize them  
# 

#%%


# insert path for db_access, a library comprising methods for working with the db
import sys
#sys.path.append('../../../btemp-db/')
import pandas as pd
import OpenWF.db_access as db_access


#%%
# relative path to the .db file, which is the actual database
db_path = '../../../../ETHeatflow/dbase_model_btemps.db'
#db_path = '../../../borehole_temp_analysis/data/interim/data_nagra_12-61/d_nab_database/dbase_model_btemps.db'


#%%
# connect to the database and get information about stored tables
conn, c = db_access.connect(db_path)


# At this point, we successfully connected to the database. One next step would be to see, what different tables are stored in the database. `db_access` provides you with methods to do so. Of course, one can directly use an SQL query to do so. For user convenience, such queries are wrapped in some python methods of `db_access`. For instance, let's check the names of tables in the database:

#%%
# SQL query
c.execute("SELECT name FROM sqlite_master WHERE type='table';")
print(c.fetchall())


#%%
# db_access method
db_access.get_tables(c)


# Essentially, these two cells of code do the same thing. In the `db_access` method, the `c.execute` and `c.fetchall` commands are bundled in one method, `.get_tables()`. The result are the three tables:  
# * borehole_information_temperatures  
# * temperature_data (with one backup table, marked with extension \_bak)  
# * sample_information_petrophysics  
# 
# In its current state, `db_access` comprises very basic query methods. More specific data-queries still need to be done via the `c.execute` and `c.fetchall` chain which is extremely versatile.  
# For instance, consider out of the over 30000 data entries, we want to get all temperature measurements for Borehole Nr. 111. 

#%%
c.execute("SELECT * FROM {tn} WHERE {idf}=111;".format(tn='temperature_data', idf='Nr'))
print(c.fetchall())


# To get the name of this borehole, we can relate to the table *borehole_information_temperatures* and query the name for the borehole with Nr. 111 in the exact same way:

#%%


c.execute("SELECT {param} FROM {tn} WHERE {idf}=111;".format(param='Namenach',
                                                             tn='borehole_information_temperatures', idf='Nr'))
print(c.fetchall())


# To know which columns are available to choose from as `{param}` in the `execute` command, we can either list names fetched by an `execute` command:

#%%


nam = c.execute("select * from borehole_information_temperatures")
names = list(map(lambda x: x[0], nam.description))
print(names)


# ... or use a `db_access` method which returns this list of table headers:

#%%


db_access.get_columns(c,table='borehole_information_temperatures')


# Now back to the query above, where we asked the database to provide all data for borehole *Riehen-1*, i.e. borehole Nr. 111. The query returns a list of table rows fitting the query command. While usable, it is difficult to read, at least for humans. This is, where pandas comes into play. As an extensive data analysis library, [pandas](https://pandas.pydata.org/) provides a lot of tools to deal with a database and present them in [dataframes](https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.html), which can be displayed in a way more organized way. Below, we submit a query for the temperature data for borehole Nr. 111 and display it.

#%%


# query database for Borehole Nr. 111 and store it in the dataframe df.
df = pd.read_sql_query("select * from temperature_data where Nr = 111;", conn)
df.head()


# Next to readability, another advantage of querying via pandas, and storing the result in a dataframe, is visualization. Pandas features some plotting functions, which can quickly plot parameters in a dataframe. For example, let's plot `Depth` versus `Temperature`:

#%%


df.plot.scatter(x='Temperature', y='Depth_asl', s=50);

#%%


c.close()
conn.close()

