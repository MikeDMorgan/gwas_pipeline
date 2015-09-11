from CGATReport.Tracker import *
import pandas as pd
from pandas.io import sql

class GenderPlotter(TrackerSQL):

    pattern = "(.+)"

    def __call__(self, track, slice=None):

        column = "f_31_0_0"

        statement = "SELECT f_eid, %(column)s from ukb4882"

        df = sql.read_sql(statement,
                          dbh,
                          index_col="f_eid")
        return df
