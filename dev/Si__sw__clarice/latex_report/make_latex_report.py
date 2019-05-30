from latex_report import LatexReport
from latex_report import insert_figure_str

class PyposmatLatexReport(LatexReport):
    def __init__(self,data_directory,iterations):
        LatexReport.__init__(self)
        self.data_directory = data_directory
        self.iterations=iterations
    def write(self,filename):
        assert isinstance(filename,str)

        with open(filename,'w') as f:
            f.write(self.get_header_section_str())
            f.write(self.get_title_str())

            f.write(self.qoi_analysis_str(
                data_directory=self.data_directory,
                iterations=self.iterations))

            if self.bibtex_fn is not None:
                f.write(self.get_bibiography_str())
            f.write(self.get_footer_str())

    def qoi_analysis_str(self,
            data_directory,
            iterations):
        s = "\\section{Qoi Analysis}\n"

        qoi_1d_plot_dir = 'qoi_1d_plots'
        s += self.qoi_1d_analysis(
                data_directory=data_directory,
                plot_directory=qoi_1d_plot_dir,
                iterations=iterations)
        return s

    def qoi_1d_analysis(self,
            data_directory,
            plot_directory,
            iterations):
        from make_qoi_plots import make_qoi_plots
        
        s = "\\subsection{1d analysis}\n"
        plot_fns = make_qoi_plots(
                data_directory=data_directory,
                plot_directory=plot_directory,
                iterations=iterations
                )
        for plot_fn in plot_fns:
            caption_s = "\\detokenize{"+plot_fn+"}"
            s += insert_figure_str(filename=plot_fn,caption=caption_s)
        return s

if __name__ == "__main__":
    import os
    import pypospack.utils

    pypospack_root_dir = pypospack.utils.get_pypospack_root_directory()
    data_dir = os.path.join(pypospack_root_dir,'data','Si__sw__data','pareto_optimization_p_3.5_q_0.5')
    iterations=[0,1,4,9,19]

    o = PyposmatLatexReport(
            data_directory=data_dir,
            iterations=iterations
            )
    o.title = "Pyposmat Simulation Report"
    o.write(filename='latex_report_test.tex')
