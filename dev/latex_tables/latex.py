class LatexTableError(Exception): pass

class LatexTable():
    DEFAULT_COLUMN_FORMAT = 'c'

    def get_column_format_string(n_columns=None,column_format_list=None):
        
        assert isinstance(n_columns,int) or isinstance(column_format_list,list)
        
        if isinstance(n_columns,int):
            col_fmt_str = n_columns*LatexTable.DEFAULT_COLUMN_FORMAT

        if isinstance(column_format_list,list):
            col_fmt_str = "".join(column_format_list)
              
        return col_fmt_str

    def get_table_string(column_formats,table_loc="ht"):

        assert isinstance(column_format,str)
        s = (
                "\\begin{table}["+table_loc+"]\n"
                "\\begin{tabular}{"+column_formats+"}\n"
                "\\hline\n"
                
                "\\hline\n"
                "\\hline\n"
                "\\end{tabular}\n"
                "\\end{table}\n"
            )
        return s

if __name__ == "__main__":
    n_columns = 5
    column_formats = LatexTable.get_column_format_string(n_columns=n_columns)
    s = LatexTable.get_table_string(column_formats=column_formats)
    print(s)
