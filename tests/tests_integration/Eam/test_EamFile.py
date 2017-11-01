import pytest

class Test__get_eam_radii_vector():
    pass

class Test__EamFile__init(object):

    def setup__variables(self):
        pass

    def test__import(self):
        import pypospack.eamtools as eamtools

    def setup__init(self):
        import pypospack.eamtools as eamtools
        self.eamfile = eamtools.EamSetflFile()

    def test__init(self):
        self.setup__init()
        # an eam file is allowed to have three comment lines
        assert isinstance(self.eamfile.comments, list)
        assert len(self.eamfile.comments) == 3

        assert self.eamfile.title is None
        assert self.eamfile.symbols is None
        assert self.eamfile.filename is None

        assert self.eamfile.r is None
        assert self.eamfile.embed is None
        assert self.eamfile.dens is None

        assert self.eamfile.N_rho is None
        assert self.eamfile.N_r == 500
        assert self.eamfile.d_r == 500

        assert self.eamfile.rcut_g == None

    def test__property_comments(self):
        self.setup__init()
        comments = ["comment line 1",
                    "comment line 2",
                    "comment line 3"]
        self.eamfile.comments = comments
        n_comments = len(comments)
        for i in range(n_comments):
            assert self.eamfile.comments[i] == comments[i]

    def test__property_comments__too_many_lines(self):
        self.setup__init()
        comments = ["comment line 1",
                    "comment line 2",
                    "comment line 3",
                    "comment line 4"]

        with pytest.raises(
                ValueError,
                message="comments must be list of length 3"):
            self.eamfile.comments = comments

    def test__property_comments__too_few_lines(self):
        self.setup__init()
        comments = ["comment line 1",
                    "comment line 2"]

        with pytest.raises(
                ValueError,
                message="comments must be list of length 3"):
            self.eamfile.comments = comments
    
    def test__property_comments__not_a_string(self):
        self.setup__init()
        comments = [1,2,3]
        self.eamfile.comments = comments
       
        n_comments = len(comments)
        for i in range(n_comments):
            assert self.eamfile.comments[i] == str(comments[i])

    def test__get_comments_to_string__no_args(self):
        self.setup__init()
        comments = ["comment 1","comment 2","comment 3"]
        comment_str = "\n".join(comments)

        self.eamfile.comments = comments
        assert self.eamfile.get_comments_as_string() == comment_str

    def test__n_symbols_line_as_string__one_element(self):
        self.setup__init()
        symbols = ['Ni']
        n_symbols_line_str = " ".join(
                [str(len(symbols))] + symbols)

        self.eamfile.symbols = symbols
        assert self.eamfile.get_n_symbols_line_as_string() \
                == n_symbols_line_str

    def test__n_symbols_line_as_string__two_elements(self):
        self.setup__init()
        symbols = ['Ni','Al']
        n_symbols_line_str = " ".join(
                [str(len(symbols))] + symbols)

        self.eamfile.symbols = symbols
        assert self.eamfile.get_n_symbols_line_as_string() \
                == n_symbols_line_str


    def test__get_eam_head(self):
        pass
    #def setup__get_eam_head(self):
    #    self.setup__init()
    #
    #def test__get_eam_head(self):
    #    self.setup__write_eam_head()
    #    self.eamfile.write_eam_head()
    #
    #def test__write_setfl(self):
    #    pass
