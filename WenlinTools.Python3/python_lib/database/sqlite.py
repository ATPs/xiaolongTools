#lib for handling database

import sqlite3 as lite
import os, sys

class sqlite(object):
    def __init__(self, f_db, tableName, quiet=True):
        self.f_db = f_db
        self.tableName = tableName
        self.quiet = quiet

    def createTable(self, columnTags):
        """
        column is a enumerable object
        example:
            createTable(column=["date text", "trans text", "symbol text", "qty real", "price real"])
        """
        cmd = "CREATE TABLE %s (%s)" % (self.tableName, ','.join(columnTags))
        self.execute(cmd, change=True)


    @property
    def table(self):
        return self.readAll()


    def readAll(self, limit=9999):
        """
        read the whole database into a list of tuple
        when reach the line limit, stop and print an warning
        """
        if not os.path.exists(self.f_db):
            return []


        conn = lite.connect(self.f_db)
        c = conn.cursor()
        cmd = "SELECT * FROM %s" % self.tableName
        c.execute(cmd)

        final = [[i[0] for i in c.description]]
        ele = True
        count = 0
        while (count < limit):
            ele = c.fetchone()
            if ele:
                final.append(ele)
                count += 1
            else:
                break
        conn.close()

        if count == limit:
            print("Warning: reach the limit! only take the first %s lines!" % limit, file=sys.stderr)
        return final


    def execute(self, cmd, fetch=False, change=False):
        """
        a simple wrapper to execute certain database command
        if fetch == True, return fectchall()
        """
        if not self.quiet:
            print("executing database command %s" % cmd, file=sys.stderr)
        conn = lite.connect(self.f_db)
        c = conn.cursor()
        c.execute(cmd)
        if change:
            conn.commit()
        if fetch:
            a = c.fetchall()
        conn.close()
        if fetch:
            return a

    def connect(self):
        """
        make the connection object
        """
        conn = lite.connect(self.f_db)
        return conn

    def insert(self, alistoflist, column=None):
        """
        if you want to insert specific column, use column=(tag_names,)
        example:
            insert([(wenlin, love), (fang, forever)])
            or
            insert([wenlin, love, Fang])
        """
        if type(alistoflist[0]) is not list and type(alistoflist[0]) is not tuple:
            alistoflist = [alistoflist]

        sample = alistoflist[0]
        values = ["?"] * len(sample)
        if column:
            cmd = "INSERT INTO %s (%s) VALUES (%s)" % (self.tableName, ",".join(column), ",".join(values))
        else:
            cmd = "INSERT INTO %s VALUES (%s)" % (self.tableName, ",".join(values))

        conn = lite.connect(self.f_db)
        c = conn.cursor()
        c.executemany(cmd, [tuple(i) for i in alistoflist])
        conn.commit()
        conn.close()

        if not self.quiet:
            print("successfully insert your record!", file=sys.stderr)


    def select(self, argv=None, column="*"):
        """
        example:
            select("where name='wenlin'", column="name, my")
        """
        cmd = "SELECT %s FROM %s" % (column, self.tableName)

        if argv:
            cmd += " %s" % argv

        a = self.execute(cmd, fetch=True, change=False)
        return a


    def delete(self, argv):
        """
        example:
            delete("where name='wenlin'")
        """
        cmd = "DELETE FROM %s %s" % (self.tableName, argv)
        self.execute(cmd, change=True)


    def update(self, where, new):
        """
        example: find the record with name=wenlin, then change it to name=Fang
        example:
            update(where="name='Wenlin'", new="name=wenlin")
        """
        cmd = "UPDATE %s SET %s WHERE %s;" % (self.tableName, new, where)
        self.execute(cmd, change=True)


    def alter(self, argv):
        """
        used to do ALTER command such as add, drop column.
        example:
            alert("ADD column_name datatype")
        """
        cmd = "ALTER TABLE %s %s;" % (self.tableName, argv)
        self.execute(cmd, change=True)






