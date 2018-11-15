#compared to mysql, basically reduced the connection times

import sys
import os#used by self.backup
import MySQLdb as mdb
import _mysql_exceptions

class mysql(object):
    """
    old setting for mysql database in pixe
    user="wenlin",
    passwd="liwenlin",
    host="pixe.local",
    unix_socket="/usr5/local/mysql/log/mysql.sock",
    port=43000,
    """
    def __init__(self,
            user="wenlin",
            passwd="FangZhang1988",
            #host="feronia.local",
            host="10.0.0.63",
            unix_socket="/data/wenlin_database/mysql/log/mysql.sock",
            port=44000,
            db="wenlindb",
            tableName="test"):
        self.user = user
        self.passwd = passwd
        self.host = host
        self.unix_socket = unix_socket
        self.port = port
        self.db = db
        self.tableName = tableName
        self.quiet = True
        self.conn = self._connect()


    def _connect(self):
        con = mdb.connect(
                user=self.user,
                passwd=self.passwd,
                host=self.host,
                unix_socket=self.unix_socket,
                port=self.port,
                db=self.db)
        return con

    def close(self):
        self.conn.close()


    def execute(self, cmd, fetch=True):
        if not self.quiet:
            print(cmd)
        con = self.conn
        cur = con.cursor()
        cur.execute(cmd)
        if fetch:
            alist = cur.fetchall()
            return alist


    def createTable(self, column=[], name="", new=False):
        """
        example: column=["first TINYTEXT, last TINYTEXT"]
        or column=["first TINYTEXT", "last TINYTEXT"]
        """
        if self.tableName != None:
            name = self.tableName
        elif name == "":
            raise IOError("Table name is not given!")
        cmd = "CREATE TABLE %s (%s);" % (name,",".join(column))
        if new == True:
            cmd = "DROP TABLE IF EXISTS %s; %s" % (name, cmd)

        try:
            self.execute(cmd)
        except _mysql_exceptions.OperationalError:
            if not self.quiet:
                print(sys.stderr, "table %s exists!" % name)


    @property
    def header(self):
        cmd = "SELECT * FROM %s LIMIT 1;" % (self.tableName)
        c = self.conn.cursor()
        c.execute(cmd)

        final = [i[0] for i in c.description]
        return final


    @property
    def table(self):
        return self.readAll()


    def readAll(self, name='', limit=9999):
        """
        read the whole database into a list of tuple
        when reach the line limit, stop and print an warning
        """
        if name == "" and self.tableName != None:
            name = self.tableName

        if self.tableName == None and name == '':
            print("Error: No table name specified!", file=sys.stderr)
            return []

        cmd = "SELECT * FROM %s LIMIT %s;" % (name, limit)
        con = self.conn
        c = con.cursor()
        c.execute(cmd)

        final = [[i[0] for i in c.description]]
        final += c.fetchall()

        if len(final) >= limit:
            print("Warning: reach the limit! only take the first %s lines!" % limit, file=sys.stderr)
        return final


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

        values = ["('%s')" % "','".join([str(i) for i in alist]) for alist in alistoflist]
        if column:
            cmd = "INSERT INTO %s (%s) VALUES %s;" % (self.tableName, ",".join(column), ",".join(values))
        else:
            cmd = "INSERT INTO %s VALUES %s;" % (self.tableName, ",".join(values))

        self.execute(cmd)

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

        a = self.execute(cmd, fetch=True)
        return a


    def delete(self, argv):
        """
        example:
            delete("where name='wenlin'")
        """
        cmd = "DELETE FROM %s %s" % (self.tableName, argv)
        self.execute(cmd)


    def update(self, cmd):
        """
        example: find the record with name=wenlin, then change it to name=Fang
        example:
            update( "SET name='Wenlin' WHERE name='wenlin'")
        """
        cmd = "UPDATE %s %s;" % (self.tableName, cmd)
        self.execute(cmd)


    def alter(self, argv):
        """
        used to do ALTER command such as add, drop column.
        example:
            alert("ADD column_name datatype")
        """
        cmd = "ALTER TABLE %s %s;" % (self.tableName, argv)
        self.execute(cmd)


    def debug(self):
        self.quiet = False

    def backup(self, dn=None, allTable=False):
        #cmd = "/data/wenlin_database/mysql/bin/mysqldump -u wenlin -p liwenlin -h feronia.local --socket=/data/wenlin_database/mysql/log/mysql.sock -P 44000 %s %s > %s" % (self.db, self.tableName, dn)
        if not dn:
            if allTable:
                dn = "%s.sql" % self.db
            else:
                dn = "%s.%s.sql" % (self.db, self.tableName)

        if allTable:
            cmd = "/data/wenlin_database/mysql/bin/mysqldump -u wenlin -p -h feronia.local --socket=/data/wenlin_database/mysql/log/mysql.sock -P 44000 %s > %s" % (self.db, dn)
        else:
            cmd = "/data/wenlin_database/mysql/bin/mysqldump -u wenlin -p -h feronia.local --socket=/data/wenlin_database/mysql/log/mysql.sock -P 44000 %s %s > %s" % (self.db, self.tableName, dn)
        os.system(cmd)
        if not self.quiet:
            if allTable:
                print("successfully backup all table in %s into %s!" % (self.db, dn))
            else:
                print("successfully backup %s.%s into %s!" % (self.db, self.tableName, dn))


class NCBI_db(mysql):
    def __init__(self,tableName):
        mysql.__init__(self)
        self.tableName = tableName
        self.db = "NCBI_refs"
        self.close()
        self.conn = self._connect()


class uniprot(mysql):
    def __init__(self, tableName):
        mysql.__init__(self, tableName=tableName, db="uniprot")




class Entrez(mysql):
    def __init__(self, tableName):
        mysql.__init__(self, tableName=tableName, db="Entrez")





class pixe_server(mysql):
    """
    old setting for mysql database in pixe
    user="wenlin",
    passwd="liwenlin",
    host="pixe.local",
    unix_socket="/usr5/local/mysql/log/mysql.sock",
    port=43000,
    """
    def __init__(self):
        mysql.__init__(self,
            user="wenlin",
            passwd="liwenlin",
            host="pixe.local",
            unix_socket="/usr5/local/mysql/log/mysql.sock",
            port=43000,
            db="test",
            tableName="test")
