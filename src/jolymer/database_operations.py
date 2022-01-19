#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct 18 15:07:19 2020

@author: johannes
"""

import sqlite3
import getpass
import pandas as pd


class dbopen:
    """
    Simple CM for sqlite3 databases. Commits everything at exit.
    """
    # path = '$HOME/LRZ Sync+Share/master-thesis/database/database.db'
    path = f"C:\\Users\\{getpass.getuser()}\\LRZ Sync+Share\\master-thesis\\database\\database.db"
    if getpass.getuser() == 'johannes':
        path = '/home/johannes/LRZ Sync+Share/master-thesis/database/database.db'
    # path = '../../database/database.db'
    print(path)

    def __init__(self):
        pass

    def __enter__(self):
        self.conn = sqlite3.connect(self.path)
        self.cursor = self.conn.cursor()
        return self.cursor

    def __exit__(self, exc_class, exc, traceback):
        self.conn.commit()
        self.conn.close()


def execute(query):
    with dbopen() as c:
        c.execute(query)
        out = c.fetchall()
    return out

def insert_values(table, args):
    qs = "?"
    for arg in range(len(args) -1):
        qs += ", ?"
    query = f"INSERT INTO {table} values ({qs});"
    with dbopen() as c:
        c.execute(query, args)
    return get_table(table)

def get_table(tablename, columns = '*'):
    query = f'SELECT {columns} FROM {tablename};'
    with dbopen() as c:
        conn = c.connection
        df = pd.read_sql(query, conn)
    return df

def print_table(*args, **kwargs):
    df = get_table(*args, **kwargs)
    print(df)

def delete_row(table, key, value):
    "deletes one row from table and returns the table ad pandas df"
    query = f'DELETE FROM {table}  WHERE {key} = {value};'
    with dbopen() as c:
         c.execute(query)
    return get_table(table)

def delete_id(table, id):
    return delete_row(table, 'id', id)

def create_generic_measurement_table(tablename):

    query = f"""CREATE TABLE {tablename}(
        id INTEGER PRIMARY KEY,
        date TEXT,
        sample TEXT,
        comment TEXT);
        """
    with dbopen() as c:
        c.execute(query)
