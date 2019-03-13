#! /usr/bin/env python
"""
Created on Fri Dec  2 14:37:17 2016

@author: hamedkhakzad
"""

# Simple code to convert CSV to SQLITE3 using pandas

def csv2db (input1, num):

    import pandas as pd
    import sqlite3
    
    table_name = 'hrMS1_result'+str(num)+'.db'

    con = sqlite3.connect(table_name)
    cur = con.cursor()
    cur.execute(
    "CREATE TABLE IF NOT EXISTS 'result' (lightmonomz REAL,heavymonomz REAL,rtapexlight REAL,rtapexheavy REAL,masslight REAL,massheavy REAL,intlight REAL,intheavy REAL,intscore REAL,rtscore REAL,massscore REAL,scoresum REAL,z REAL,xlink TEXT,model_tag INTEGER,iteration INTEGER,ratio INTEGER);")
    
    result = pd.DataFrame(input1)
    result.to_sql("result", con, if_exists='append', index=False)
    
    # Insert a row of data
    cur.execute("SELECT DISTINCT xlink, COUNT(xlink) AS cnt FROM result GROUP BY xlink")
    rows = cur.fetchall()
    # update `iteration` value
    cur.executemany("UPDATE result SET iteration=? WHERE xlink=?",
                  [(row[1], row[0]) for row in rows])
    
    cur.execute("SELECT DISTINCT xlink, count(case when model_tag='1' then 1 else null end) as ratio from result GROUP BY xlink")
    rows2 = cur.fetchall()
    # update `iteration` value
    cur.executemany("UPDATE result SET ratio=? WHERE xlink=?",
                  [(row[1], row[0]) for row in rows2])
    # Save (commit) the changes
    con.commit()
    
    # Extract top XLs
    cur.execute("SELECT DISTINCT xlink from result where model_tag='1'")
    top_XL = cur.fetchall()
    
    # We can also close the connection if we are done with it.
    # Just be sure any changes have been committed or they will be lost.
    con.close()
    
    return top_XL
