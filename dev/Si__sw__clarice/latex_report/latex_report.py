from pathlib import Path

def insert_figure_str(filename,caption=None):

    s = "\\begin{figure}\n"
    s += "  \\includegraphics{"+str(Path(filename).with_suffix(''))+"}\n"
    s += "  \\centering\n"
    if caption is not None:
        s += "  \\caption{"+caption+"}\n"
    s += "\\end{figure}\n"
    return s

class LatexReport(object):

    def __init__(self):

        self.document_type='article'
        self.packages=[
                'amsmath',
                'cite',
                'graphicx']
        self.title = 'no title'
        self.author = 'generated by pyposmat'
        self.date = '\\today'

        self.bibtex_fn = None
        self.bibtex_style = 'acm'
    
    def write(self,filename):
        assert isinstance(filename,str)

        with open(filename,'w') as f:
            f.write(self.get_header_section_str())
            f.write(self.get_title_str())
            if self.bibtex_fn is not None:
                f.write(self.get_bibiography_str())
            f.write(self.get_footer_str())

    def get_header_section_str(self):
        s = "\\documentclass{"+self.document_type+"}\n"
        for package in self.packages:
            s += "\\usepackage{"+package+"}\n"
            s += "\\usepackage[T1]{fontenc}"
        s += "\\begin{document}\n"
        return s

    def get_title_str(self):
        s = "\\title{"+self.title+"}\n"
        s += "\\author{"+self.author+"}\n"
        s += "\\date{"+self.date+"}\n"
        s += "\\maketitle\n"
        return s

    def get_bibliography_str(self):
        s = "\\bibliography{"+self.bibtex_fn+"}\n"
        s += "\\bibliographystyle{"+self.bibtex_style+"}\n"
        return s

    def get_footer_str(self):
        s = "\\end{document}"
        return s
if __name__ == "__main__":
    o = LatexReport()
    o.write(filename='latex_report_test.tex')
