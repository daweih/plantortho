Bio::Phylo::IO(3)     User Contributed Perl Documentation    Bio::Phylo::IO(3)

NAME
       Bio::Phylo::IO - Front end for parsers and serializers

SYNOPSIS
        use Bio::Phylo::IO qw(parse unparse);

        # returns an unblessed array reference of block objects,
        # i.e. taxa, matrix or forest objects
        my $blocks = parse(
           '-file'     => $file,
           '-format'   => 'nexus',
           '-encoding' => ':encoding(UTF-8)', # optional, default is system-dependent
        );

        for my $block ( @{ $blocks } ) {
           if ( $block->isa('Bio::Phylo::Taxa') ) {
               my $taxa = $block;
               # do something with the taxa
           }
        }

        # returns a Bio::Phylo::Project object
        my $project = parse(
               '-file'       => $file,
               '-format'     => 'nexus',
               '-as_project' => 1
        )
        my ($taxa) = @{ $project->get_taxa };

        # parsing a tree from a newick string
        my $tree_string = '(((A,B),C),D);';
        my $tree = Bio::Phylo::IO->parse(
           '-string' => $tree_string,
           '-format' => 'newick',
        )->first;

        # note: newick parsers return
        # 'Bio::Phylo::Forest'! Call
        # ->first to retrieve the first
        # tree of the forest.

        # prints 'Bio::Phylo::Forest::Tree'
        print ref $tree, "\n";

        # if the tree is very large and you need only some terminal nodes from it
        $simplified_tree = Bio::Phylo::IO->parse(
           '-string' => $tree_string,
           '-format' => 'newick',
           '-keep'   => ['A', 'D'], # nodes to keep
           '-ignore_comments' => 1, # treats [] symbols as part of taxon name
        )->first;

        # parsing a table
        my $table_string = qq(A,1,2|B,1,2|C,2,2|D,2,1);
        my $matrix = Bio::Phylo::IO->parse(
           '-string'   => $table_string,
           '-format'   => 'table',

           # Data type, see Bio::Phylo::Parsers::Table
           '-type'     => 'STANDARD',

           # field separator
           '-fieldsep' => ',',

           # line separator
           '-linesep'  => '|'
        );

        # prints 'Bio::Phylo::Matrices::Matrix'
        print ref $matrix, "\n";

        # parsing a list of taxa
        my $taxa_string = 'A:B:C:D';
        my $taxa = Bio::Phylo::IO->parse(
           '-string'   => $taxa_string,
           '-format'   => 'taxlist',
           '-fieldsep' => ':'
        );

        # prints 'Bio::Phylo::Taxa'
        print ref $taxa, "\n";

        # matches taxon names in tree to $taxa object
        $tree->cross_reference($taxa);

        # likewise for matrix
        $matrix->cross_reference($taxa);

        print unparse(

           # pass the tree object,
           # crossreferenced to taxa, which
           # are crossreferenced to the matrix
           '-phylo' => $tree,
           '-format' => 'pagel'
        );

        # prints a pagel data file:
        #4 2
        #A,n1,0.000000,1,2
        #B,n1,0.000000,1,2
        #n1,n2,0.000000
        #C,n2,0.000000,2,2
        #n2,n3,0.000000
        #D,n3,0.000000,2,1

DESCRIPTION
       The IO module is the front end for parsing and serializing phylogenetic
       data objects. It is a non-OO module that optionally exports the 'parse'
       and 'unparse' subroutines into the caller's namespace, using the "use
       Bio::Phylo::IO qw(parse unparse);" directive. Alternatively, you can
       call the subroutines as class methods. The "parse" and "unparse"
       subroutines load and dispatch the appropriate sub-modules at runtime,
       depending on the '-format' argument.

   CLASS METHODS
       parse()
           Parses a file or string.

            Type    : Class method
            Title   : parse
            Usage   : my $obj = Bio::Phylo::IO->parse(%options);
            Function: Creates (file) handle,
                      instantiates appropriate parser.
            Returns : A Bio::Phylo::* object
            Args    : -file    => (path),
                       or
                      -string  => (scalar),
                      or
                      -handle  => (IO::Handle object)
                      or
                      -url     => (url string)
                      -format  => (description format),
                      -(other) => (parser specific options)
            Comments: The parse method makes assumptions about
                      the capabilities of Bio::Phylo::Parsers::*
                      modules: i) their names match those of the
                      -format => (blah) arguments, insofar that
                      ucfirst(blah) . '.pm' is an existing module;
                      ii) the modules implement a _from_handle,
                      or a _from_string method. Exceptions are
                      thrown if either assumption is violated.

                      If @ARGV contains even key/value pairs such
                      as "format newick file <filename>" (note: no
                      dashes) these will be prepended to @_, for
                      one-liners.

       parse_matrix()
           Parses a file or string.

            Type    : Class method
            Title   : parse_matrix
            Usage   : my $matrix = Bio::Phylo::IO->parse_matrix(%options);
            Function: Creates (file) handle,
                      instantiates appropriate parser.
            Returns : A Bio::Phylo::Matrices::Matrix object
            Args    : Same as parse()
            Comments: This method is syntactical sugar to get the first matrix
                      out of a file/handle/string

       parse_tree()
           Parses a file or string.

            Type    : Class method
            Title   : parse_tree
            Usage   : my $tree = Bio::Phylo::IO->parse_tree(%options);
            Function: Creates (file) handle,
                      instantiates appropriate parser.
            Returns : A Bio::Phylo::Forest::Tree object
            Args    : Same as parse()
            Comments: This method is syntactical sugar to get the first tree
                      out of a file/handle/string

       unparse()
           Unparses object(s) to a string.

            Type    : Class method
            Title   : unparse
            Usage   : my $string = Bio::Phylo::IO->unparse(
                          %options
                      );
            Function: Turns Bio::Phylo object into a
                      string according to specified format. If an
                      optional -file or -handle argument is provided
                      the string is also written to that.
            Returns : SCALAR
            Args    : -phylo   => (Bio::Phylo object),
                      -format  => (description format),
                      -(other) => (parser specific options)
                      -file    => (optional: a file path to open and write to)
                      or
                      -handle  => (optional: a handle to write to)

       can_read()
           Tests whether Bio::Phylo::IO can read provided syntax format.

            Type    : Class method
            Title   : can_read
            Usage   : &do_something if Bio::Phylo::IO->can_read('foo');
            Function: Tests whether Bio::Phylo::IO can read provided syntax format.
            Returns : Boolean
            Args    : A syntax format name, like "nexml"

       can_write()
           Tests whether Bio::Phylo::IO can write provided syntax format.

            Type    : Class method
            Title   : can_write
            Usage   : &do_something if Bio::Phylo::IO->can_write('foo');
            Function: Tests whether Bio::Phylo::IO can write provided syntax format.
            Returns : Boolean
            Args    : A syntax format name, like "nexml"

SEE ALSO
       There is a mailing list at
       https://groups.google.com/forum/#!forum/bio-phylo
       <https://groups.google.com/forum/#!forum/bio-phylo> for any user or
       developer questions and discussions.

       Bio::Phylo::Parsers::Fasta
       Bio::Phylo::Parsers::Newick
       Bio::Phylo::Parsers::Nexml
       Bio::Phylo::Parsers::Nexus
       Bio::Phylo::Parsers::Phylip
       Bio::Phylo::Parsers::Phyloxml
       Bio::Phylo::Parsers::Table
       Bio::Phylo::Parsers::Taxlist
       Bio::Phylo::Parsers::Tolweb
       Bio::Phylo::Unparsers::Mrp
       Bio::Phylo::Unparsers::Newick
       Bio::Phylo::Unparsers::Nexml
       Bio::Phylo::Unparsers::Nexus
       Bio::Phylo::Unparsers::Pagel
       Bio::Phylo::Unparsers::Phylip
       Bio::Phylo::Unparsers::Phyloxml
       Bio::Phylo::Manual
           Also see the manual: Bio::Phylo::Manual and
           <http://rutgervos.blogspot.com>

CITATION
       If you use Bio::Phylo in published research, please cite it:

       Rutger A Vos, Jason Caravas, Klaas Hartmann, Mark A Jensen and Chase
       Miller, 2011. Bio::Phylo - phyloinformatic analysis using Perl.  BMC
       Bioinformatics 12:63.  http://dx.doi.org/10.1186/1471-2105-12-63
       <http://dx.doi.org/10.1186/1471-2105-12-63>

perl v5.16.3                      2014-02-08                 Bio::Phylo::IO(3)
