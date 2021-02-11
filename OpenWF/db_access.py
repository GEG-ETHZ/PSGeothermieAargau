"""db_access is a compendium of high-level python functions to mess around with a database."""

import sqlite3
import pandas as pd

def connect(db_file):
    """
    Connect to an existing sqlite database file.

    :param db_file: string - path to database file
    :return: connection object or None
    """
    try:
        conn = sqlite3.connect(db_file)
        c = conn.cursor()
        return conn, c
    except sqlite3.Error as err:
        print(err)

def close(conn, commit='Yes'):
    """Close connection to the database, commit or dont."""
    if commit == 'Yes':
        conn.commit()
    conn.close()

def get_tables(c, verbose=False):
    """Get all table names in the database."""
    c.execute("SELECT name FROM sqlite_master WHERE type='table';")
    tabs = c.fetchall()
    if verbose:
        print(tabs)
    return tabs

def get_columns(c, table, verbose=False):
    """Get all columns in a specified table."""
    head = c.execute("select * from " + table)
    names = list(map(lambda x: x[0], head.description))
    if verbose:
        print(names)
    return(names)
