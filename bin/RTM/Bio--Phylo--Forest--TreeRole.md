Bio::Phylo::Forest::TrUserlContributed Perl DocBio::Phylo::Forest::TreeRole(3)

NAME
       Bio::Phylo::Forest::TreeRole - Extra behaviours for a phylogenetic tree

SYNOPSIS
        # some way to get a tree
        use Bio::Phylo::IO;
        my $string = '((A,B),C);';
        my $forest = Bio::Phylo::IO->parse(
           -format => 'newick',
           -string => $string
        );
        my $tree = $forest->first;

        # do something:
        print $tree->calc_imbalance;

        # prints "1"

DESCRIPTION
       The object models a phylogenetic tree, a container of
       Bio::Phylo::Forest::Node objects. The tree object inherits from
       Bio::Phylo::Listable, so look there for more methods.

METHODS
   CONSTRUCTORS
       new()
           Tree constructor.

            Type    : Constructor
            Title   : new
            Usage   : my $tree = Bio::Phylo::Forest::Tree->new;
            Function: Instantiates a Bio::Phylo::Forest::Tree object.
            Returns : A Bio::Phylo::Forest::Tree object.
            Args    : No required arguments.

       new_from_bioperl()
           Tree constructor from Bio::Tree::TreeI argument.

            Type    : Constructor
            Title   : new_from_bioperl
            Usage   : my $tree =
                      Bio::Phylo::Forest::Tree->new_from_bioperl(
                          $bptree
                      );
            Function: Instantiates a
                      Bio::Phylo::Forest::Tree object.
            Returns : A Bio::Phylo::Forest::Tree object.
            Args    : A tree that implements Bio::Tree::TreeI

   QUERIES
       get_midpoint()
           Gets node that divides tree into two distance-balanced partitions.

            Type    : Query
            Title   : get_midpoint
            Usage   : my $midpoint = $tree->get_midpoint;
            Function: Gets node nearest to the middle of the longest path
            Returns : A Bio::Phylo::Forest::Node object.
            Args    : NONE
            Comments: This algorithm was ported from ETE

       get_terminals()
           Get terminal nodes.

            Type    : Query
            Title   : get_terminals
            Usage   : my @terminals = @{ $tree->get_terminals };
            Function: Retrieves all terminal nodes in
                      the Bio::Phylo::Forest::Tree object.
            Returns : An array reference of
                      Bio::Phylo::Forest::Node objects.
            Args    : NONE
            Comments: If the tree is valid, this method
                      retrieves the same set of nodes as
                      $node->get_terminals($root). However,
                      because there is no recursion it may
                      be faster. Also, the node method by
                      the same name does not see orphans.

       get_internals()
           Get internal nodes.

            Type    : Query
            Title   : get_internals
            Usage   : my @internals = @{ $tree->get_internals };
            Function: Retrieves all internal nodes
                      in the Bio::Phylo::Forest::Tree object.
            Returns : An array reference of
                      Bio::Phylo::Forest::Node objects.
            Args    : NONE
            Comments: If the tree is valid, this method
                      retrieves the same set of nodes as
                      $node->get_internals($root). However,
                      because there is no recursion it may
                      be faster. Also, the node method by
                      the same name does not see orphans.

       get_cherries()
           Get all cherries, i.e. nodes that have two terminal children

            Type    : Query
            Title   : get_cherries
            Usage   : my @cherries = @{ $tree->get_cherries };
            Function: Returns an array ref of cherries
            Returns : ARRAY
            Args    : NONE

       get_all_rootings()
           Gets a forest of all rooted versions of the invocant tree.

            Type    : Query
            Title   : get_all_rootings
            Usage   : my $forest = $tree->get_all_rootings;
            Function: Returns an array ref of cherries
            Returns : Bio::Phylo::Forest object
            Args    : NONE
            Comments: This method assumes the invocant tree has a basal trichotomy.
                      "Rooted" trees with a basal bifurcation will give strange
                      results.

       get_root()
           Get root node.

            Type    : Query
            Title   : get_root
            Usage   : my $root = $tree->get_root;
            Function: Returns the root node.
            Returns : Bio::Phylo::Forest::Node
            Args    : NONE

       get_ntax()
           Gets number of tips

            Type    : Query
            Title   : get_ntax
            Usage   : my $ntax = $tree->get_ntax;
            Function: Calculates the number of terminal nodes
            Returns : Int
            Args    : NONE

       get_tallest_tip()
           Retrieves the node furthest from the root.

            Type    : Query
            Title   : get_tallest_tip
            Usage   : my $tip = $tree->get_tallest_tip;
            Function: Retrieves the node furthest from the
                      root in the current Bio::Phylo::Forest::Tree
                      object.
            Returns : Bio::Phylo::Forest::Node
            Args    : NONE
            Comments: If the tree has branch lengths, the tallest tip is
                      based on root-to-tip path length, else it is based
                      on number of nodes to root

       get_nodes_for_taxa()
           Gets node objects for the supplied taxon objects

            Type    : Query
            Title   : get_nodes_for_taxa
            Usage   : my @nodes = @{ $tree->get_nodes_for_taxa(\@taxa) };
            Function: Gets node objects for the supplied taxon objects
            Returns : array ref of Bio::Phylo::Forest::Node objects
            Args    : A reference to an array of Bio::Phylo::Taxa::Taxon objects
                      or a Bio::Phylo::Taxa object

       get_mrca()
           Get most recent common ancestor of argument nodes.

            Type    : Query
            Title   : get_mrca
            Usage   : my $mrca = $tree->get_mrca(\@nodes);
            Function: Retrieves the most recent
                      common ancestor of \@nodes
            Returns : Bio::Phylo::Forest::Node
            Args    : A reference to an array of
                      Bio::Phylo::Forest::Node objects
                      in $tree.

   TESTS
       is_binary()
           Test if tree is bifurcating.

            Type    : Test
            Title   : is_binary
            Usage   : if ( $tree->is_binary ) {
                         # do something
                      }
            Function: Tests whether the invocant
                      object is bifurcating.
            Returns : BOOLEAN
            Args    : NONE

       is_ultrametric()
           Test if tree is ultrametric.

            Type    : Test
            Title   : is_ultrametric
            Usage   : if ( $tree->is_ultrametric(0.01) ) {
                         # do something
                      }
            Function: Tests whether the invocant is
                      ultrametric.
            Returns : BOOLEAN
            Args    : Optional margin between pairwise
                      comparisons (default = 0).
            Comments: The test is done by performing
                      all pairwise comparisons for
                      root-to-tip path lengths. Since many
                      programs introduce rounding errors
                      in branch lengths the optional argument is
                      available to test TRUE for nearly
                      ultrametric trees. For example, a value
                      of 0.01 indicates that no pairwise
                      comparison may differ by more than 1%.
                      Note: behaviour is undefined for
                      negative branch lengths.

       is_monophyletic()
           Tests if first argument (node array ref) is monophyletic with
           respect to second argument.

            Type    : Test
            Title   : is_monophyletic
            Usage   : if ( $tree->is_monophyletic(\@tips, $node) ) {
                         # do something
                      }
            Function: Tests whether the set of \@tips is
                      monophyletic w.r.t. $outgroup.
            Returns : BOOLEAN
            Args    : A reference to a list of nodes, and a node.
            Comments: This method is essentially the
                      same as
                      &Bio::Phylo::Forest::Node::is_outgroup_of.

       is_paraphyletic()
            Type    : Test
            Title   : is_paraphyletic
            Usage   : if ( $tree->is_paraphyletic(\@nodes,$node) ){ }
            Function: Tests whether or not a given set of nodes are paraphyletic
                      (representing the full clade) given an outgroup
            Returns : [-1,0,1] , -1 if the group is not monophyletic
                                  0 if the group is not paraphyletic
                                  1 if the group is paraphyletic
            Args    : Array ref of node objects which are in the tree,
                      Outgroup to compare the nodes to

       is_clade()
           Tests if argument (node array ref) forms a clade.

            Type    : Test
            Title   : is_clade
            Usage   : if ( $tree->is_clade(\@tips) ) {
                         # do something
                      }
            Function: Tests whether the set of
                      \@tips forms a clade
            Returns : BOOLEAN
            Args    : A reference to an array of Bio::Phylo::Forest::Node objects, or a
                      reference to an array of Bio::Phylo::Taxa::Taxon objects, or a
                      Bio::Phylo::Taxa object
            Comments:

       is_cladogram()
           Tests if tree is a cladogram (i.e. no branch lengths)

            Type    : Test
            Title   : is_cladogram
            Usage   : if ( $tree->is_cladogram() ) {
                         # do something
                      }
            Function: Tests whether the tree is a
                      cladogram (i.e. no branch lengths)
            Returns : BOOLEAN
            Args    : NONE
            Comments:

   CALCULATIONS
       calc_branch_length_distance()
           Calculates the Euclidean branch length distance between two trees.

            Type    : Calculation
            Title   : calc_branch_length_distance
            Usage   : my $distance =
                      $tree1->calc_branch_length_distance($tree2);
            Function: Calculates the Euclidean branch length distance between two trees
            Returns : SCALAR, number
            Args    : NONE

       calc_branch_length_score()
           Calculates the squared Euclidean branch length distance between two
           trees.

            Type    : Calculation
            Title   : calc_branch_length_score
            Usage   : my $score =
                      $tree1->calc_branch_length_score($tree2);
            Function: Calculates the squared Euclidean branch
                      length distance between two trees
            Returns : SCALAR, number
            Args    : NONE

       calc_tree_length()
           Calculates the sum of all branch lengths.

            Type    : Calculation
            Title   : calc_tree_length
            Usage   : my $tree_length =
                      $tree->calc_tree_length;
            Function: Calculates the sum of all branch
                      lengths (i.e. the tree length).
            Returns : FLOAT
            Args    : NONE

       calc_tree_height()
           Calculates the height of the tree.

            Type    : Calculation
            Title   : calc_tree_height
            Usage   : my $tree_height =
                      $tree->calc_tree_height;
            Function: Calculates the height
                      of the tree.
            Returns : FLOAT
            Args    : NONE
            Comments: For ultrametric trees this
                      method returns the height, but
                      this is done by averaging over
                      all root-to-tip path lengths, so
                      for additive trees the result
                      should consequently be interpreted
                      differently.

       calc_number_of_nodes()
           Calculates the number of nodes.

            Type    : Calculation
            Title   : calc_number_of_nodes
            Usage   : my $number_of_nodes =
                      $tree->calc_number_of_nodes;
            Function: Calculates the number of
                      nodes (internals AND terminals).
            Returns : INT
            Args    : NONE

       calc_number_of_terminals()
           Calculates the number of terminal nodes.

            Type    : Calculation
            Title   : calc_number_of_terminals
            Usage   : my $number_of_terminals =
                      $tree->calc_number_of_terminals;
            Function: Calculates the number
                      of terminal nodes.
            Returns : INT
            Args    : NONE

       calc_number_of_internals()
           Calculates the number of internal nodes.

            Type    : Calculation
            Title   : calc_number_of_internals
            Usage   : my $number_of_internals =
                      $tree->calc_number_of_internals;
            Function: Calculates the number
                      of internal nodes.
            Returns : INT
            Args    : NONE

       calc_number_of_cherries()
           Calculates the number of cherries, i.e. the number of nodes that
           subtend exactly two tips. See for applications of this metric:
           http://dx.doi.org/10.1016/S0025-5564(99)00060-7
           <http://dx.doi.org/10.1016/S0025-5564(99)00060-7>

            Type    : Calculation
            Title   : calc_number_of_cherries
            Usage   : my $number_of_cherries =
                      $tree->calc_number_of_cherries;
            Function: Calculates the number of cherries
            Returns : INT
            Args    : NONE

       calc_total_paths()
           Calculates the sum of all root-to-tip path lengths.

            Type    : Calculation
            Title   : calc_total_paths
            Usage   : my $total_paths =
                      $tree->calc_total_paths;
            Function: Calculates the sum of all
                      root-to-tip path lengths.
            Returns : FLOAT
            Args    : NONE

       calc_redundancy()
           Calculates the amount of shared (redundant) history on the total.

            Type    : Calculation
            Title   : calc_redundancy
            Usage   : my $redundancy =
                      $tree->calc_redundancy;
            Function: Calculates the amount of shared
                      (redundant) history on the total.
            Returns : FLOAT
            Args    : NONE
            Comments: Redundancy is calculated as
            1 / ( treelength - height / ( ntax * height - height ) )

       calc_imbalance()
           Calculates Colless' coefficient of tree imbalance.

            Type    : Calculation
            Title   : calc_imbalance
            Usage   : my $imbalance = $tree->calc_imbalance;
            Function: Calculates Colless' coefficient
                      of tree imbalance.
            Returns : FLOAT
            Args    : NONE
            Comments: As described in Colless, D.H., 1982.
                      The theory and practice of phylogenetic
                      systematics. Systematic Zoology 31(1): 100-104

       calc_i2()
           Calculates I2 imbalance.

            Type    : Calculation
            Title   : calc_i2
            Usage   : my $ci2 = $tree->calc_i2;
            Function: Calculates I2 imbalance.
            Returns : FLOAT
            Args    : NONE
            Comments:

       calc_gamma()
           Calculates the Pybus gamma statistic.

            Type    : Calculation
            Title   : calc_gamma
            Usage   : my $gamma = $tree->calc_gamma();
            Function: Calculates the Pybus gamma statistic
            Returns : FLOAT
            Args    : NONE
            Comments: As described in Pybus, O.G. and
                      Harvey, P.H., 2000. Testing
                      macro-evolutionary models using
                      incomplete molecular phylogenies.
                      Proc. R. Soc. Lond. B 267, 2267-2272

       calc_fiala_stemminess()
           Calculates stemminess measure of Fiala and Sokal (1985).

            Type    : Calculation
            Title   : calc_fiala_stemminess
            Usage   : my $fiala_stemminess =
                      $tree->calc_fiala_stemminess;
            Function: Calculates stemminess measure
                      Fiala and Sokal (1985).
            Returns : FLOAT
            Args    : NONE
            Comments: As described in Fiala, K.L. and
                      R.R. Sokal, 1985. Factors
                      determining the accuracy of
                      cladogram estimation: evaluation
                      using computer simulation.
                      Evolution, 39: 609-622

       calc_rohlf_stemminess()
           Calculates stemminess measure from Rohlf et al. (1990).

            Type    : Calculation
            Title   : calc_rohlf_stemminess
            Usage   : my $rohlf_stemminess =
                      $tree->calc_rohlf_stemminess;
            Function: Calculates stemminess measure
                      from Rohlf et al. (1990).
            Returns : FLOAT
            Args    : NONE
            Comments: As described in Rohlf, F.J.,
                      W.S. Chang, R.R. Sokal, J. Kim,
                      1990. Accuracy of estimated
                      phylogenies: effects of tree
                      topology and evolutionary model.
                      Evolution, 44(6): 1671-1684

       calc_resolution()
           Calculates tree resolution.

            Type    : Calculation
            Title   : calc_resolution
            Usage   : my $resolution =
                      $tree->calc_resolution;
            Function: Calculates the number
                      of internal nodes over the
                      total number of internal nodes
                      on a fully bifurcating
                      tree of the same size.
            Returns : FLOAT
            Args    : NONE

       calc_branching_times()
           Calculates cumulative branching times.

            Type    : Calculation
            Title   : calc_branching_times
            Usage   : my $branching_times =
                      $tree->calc_branching_times;
            Function: Returns a two-dimensional array.
                      The first dimension consists of
                      the "records", so that in the
                      second dimension $AoA[$first][0]
                      contains the internal node references,
                      and $AoA[$first][1] the branching
                      time of the internal node. The
                      records are orderered from root to
                      tips by time from the origin.
            Returns : SCALAR[][] or FALSE
            Args    : NONE

       calc_waiting_times()
           Calculates intervals between splits.

            Type    : Calculation
            Title   : calc_waiting_times
            Usage   : my $waitings =
                      $tree->calc_waiting_times;
            Function: Returns a two-dimensional array.
                      The first dimension consists of
                      the "records", so that in the
                      second dimension $AoA[$first][0]
                      contains the internal node references,
                      and $AoA[$first][1] the waiting
                      time of the internal node. The
                      records are orderered from root to
                      tips by time from the origin.
            Returns : SCALAR[][] or FALSE
            Args    : NONE

       calc_node_ages()
           Calculates node ages.

            Type    : Calculation
            Title   : calc_node_ages
            Usage   : $tree->calc_node_ages;
            Function: Calculates the age of all the nodes in the tree (i.e. the distance
                      from the tips) and assigns these to the 'age' slot, such that,
                      after calling this method, the age of any one node can be retrieved
                      by calling $node->get_generic('age');
            Returns : The invocant
            Args    : NONE
            Comments: This method computes, in a sense, the opposite of
                      calc_branching_times: here, we compute the distance from the tips
                      (i.e. how long ago the split occurred), whereas calc_branching_times
                      calculates the distance from the root.

       calc_ltt()
           Calculates lineage-through-time data points.

            Type    : Calculation
            Title   : calc_ltt
            Usage   : my $ltt = $tree->calc_ltt;
            Function: Returns a two-dimensional array.
                      The first dimension consists of the
                      "records", so that in the second
                      dimension $AoA[$first][0] contains
                      the internal node references, and
                      $AoA[$first][1] the branching time
                      of the internal node, and $AoA[$first][2]
                      the cumulative number of lineages over
                      time. The records are orderered from
                      root to tips by time from the origin.
            Returns : SCALAR[][] or FALSE
            Args    : NONE
## NOTE:
       calc_symdiff()
           Calculates the symmetric difference metric between invocant and
           argument. This metric is identical to the Robinson-Foulds tree
           comparison distance. See
           http://dx.doi.org/10.1016/0025-5564(81)90043-2
           <http://dx.doi.org/10.1016/0025-5564(81)90043-2>

            Type    : Calculation
            Title   : calc_symdiff
            Usage   : my $symdiff =
                      $tree->calc_symdiff($other_tree);
            Function: Returns the symmetric difference
                      metric between $tree and $other_tree,
                      sensu Penny and Hendy, 1985.
            Returns : SCALAR
            Args    : A Bio::Phylo::Forest::Tree object
            Comments: Trees in comparison must span
                      the same set of terminal taxa
                      or results are meaningless.

       calc_fp()
           Calculates the Fair Proportion value for each terminal.

            Type    : Calculation
            Title   : calc_fp
            Usage   : my $fp = $tree->calc_fp();
            Function: Returns the Fair Proportion
                      value for each terminal
            Returns : HASHREF
            Args    : NONE

       calc_es()
           Calculates the Equal Splits value for each terminal

            Type    : Calculation
            Title   : calc_es
            Usage   : my $es = $tree->calc_es();
            Function: Returns the Equal Splits value for each terminal
            Returns : HASHREF
            Args    : NONE

       calc_pe()
           Calculates the Pendant Edge value for each terminal.

            Type    : Calculation
            Title   : calc_pe
            Usage   : my $es = $tree->calc_pe();
            Function: Returns the Pendant Edge value for each terminal
            Returns : HASHREF
            Args    : NONE

       calc_shapley()
           Calculates the Shapley value for each terminal.

            Type    : Calculation
            Title   : calc_shapley
            Usage   : my $es = $tree->calc_shapley();
            Function: Returns the Shapley value for each terminal
            Returns : HASHREF
            Args    : NONE

   VISITOR METHODS
       The following methods are a - not entirely true-to-form -
       implementation of the Visitor design pattern: the nodes in a tree are
       visited, and rather than having an object operate on them, a set of
       code references is used. This can be used, for example, to serialize a
       tree to a string format. To create a newick string without branch
       lengths you would use something like this (there is a more powerful
       'to_newick' method, so this is just an example):

        $tree->visit_depth_first(
               '-pre_daughter'   => sub { print '('             },
               '-post_daughter'  => sub { print ')'             },
               '-in'             => sub { print shift->get_name },
               '-pre_sister'     => sub { print ','             },
        );
        print ';';

       visit_depth_first()
           Visits nodes depth first

            Type    : Visitor method
            Title   : visit_depth_first
            Usage   : $tree->visit_depth_first( -pre => sub{ ... }, -post => sub { ... } );
            Function: Visits nodes in a depth first traversal, executes subs
            Returns : $tree
             Args    : Optional handlers in the order in which they would be executed on an internal node:

                                   # first event handler, is executed when node is reached in recursion
                                   -pre            => sub { print "pre: ",            shift->get_name, "\n" },

                                   # is executed if node has a daughter, but before that daughter is processed
                                   -pre_daughter   => sub { print "pre_daughter: ",   shift->get_name, "\n" },

                                   # is executed if node has a daughter, after daughter has been processed
                                   -post_daughter  => sub { print "post_daughter: ",  shift->get_name, "\n" },

                                   # is executed whether or not node has sisters, if it does have sisters
                                   # they're processed first
                                   -in             => sub { print "in: ",             shift->get_name, "\n" },

                                   # is executed if node has a sister, before sister is processed
                                   -pre_sister     => sub { print "pre_sister: ",     shift->get_name, "\n" },

                                   # is executed if node has a sister, after sister is processed
                                   -post_sister    => sub { print "post_sister: ",    shift->get_name, "\n" },

                                   # is executed last
                                   -post           => sub { print "post: ",           shift->get_name, "\n" },

                                   # specifies traversal order, default 'ltr' means first_daugher -> next_sister
                                   # traversal, alternate value 'rtl' means last_daughter -> previous_sister traversal
                                   -order          => 'ltr', # ltr = left-to-right, 'rtl' = right-to-left
            Comments:

       visit_breadth_first()
           Visits nodes breadth first

            Type    : Visitor method
            Title   : visit_breadth_first
            Usage   : $tree->visit_breadth_first( -pre => sub{ ... }, -post => sub { ... } );
            Function: Visits nodes in a breadth first traversal, executes handlers
            Returns : $tree
            Args    : Optional handlers in the order in which they would be executed on an internal node:

                                   # first event handler, is executed when node is reached in recursion
                                   -pre            => sub { print "pre: ",            shift->get_name, "\n" },

                                   # is executed if node has a sister, before sister is processed
                                   -pre_sister     => sub { print "pre_sister: ",     shift->get_name, "\n" },

                                   # is executed if node has a sister, after sister is processed
                                   -post_sister    => sub { print "post_sister: ",    shift->get_name, "\n" },

                                   # is executed whether or not node has sisters, if it does have sisters
                                   # they're processed first
                                   -in             => sub { print "in: ",             shift->get_name, "\n" },

                                   # is executed if node has a daughter, but before that daughter is processed
                                   -pre_daughter   => sub { print "pre_daughter: ",   shift->get_name, "\n" },

                                   # is executed if node has a daughter, after daughter has been processed
                                   -post_daughter  => sub { print "post_daughter: ",  shift->get_name, "\n" },

                                   # is executed last
                                   -post           => sub { print "post: ",           shift->get_name, "\n" },

                                   # specifies traversal order, default 'ltr' means first_daugher -> next_sister
                                   # traversal, alternate value 'rtl' means last_daughter -> previous_sister traversal
                                   -order          => 'ltr', # ltr = left-to-right, 'rtl' = right-to-left
            Comments:

       visit_level_order()
           Visits nodes in a level order traversal.

            Type    : Visitor method
            Title   : visit_level_order
            Usage   : $tree->visit_level_order( sub{...} );
            Function: Visits nodes in a level order traversal, executes sub
            Returns : $tree
            Args    : A subroutine reference that operates on visited nodes.
            Comments:

   TREE MANIPULATION
       chronompl()
           Modifies branch lengths using the mean path lengths method of
           Britton et al. (2002). For more about this method, see:
           http://dx.doi.org/10.1016/S1055-7903(02)00268-3
           <http://dx.doi.org/10.1016/S1055-7903(02)00268-3>

            Type    : Tree manipulator
            Title   : chronompl
            Usage   : $tree->chronompl;
            Function: Makes tree ultrametric using MPL method
            Returns : The modified, now ultrametric invocant.
            Args    : NONE
            Comments:

       grafenbl()
           Computes and assigns branch lengths using Grafen's method, which
           makes node ages proportional to clade size. For more about this
           method, see: <http://dx.doi.org/10.1098/rstb.1989.0106>

            Type    : Tree manipulator
            Title   : grafenbl
            Usage   : $tree->grafenbl;
            Function: Assigns branch lengths using Grafen's method
            Returns : The modified, now ultrametric invocant.
            Args    : Optional, a power ('rho') to which all node ages are raised
            Comments:

       agetobl()
           Converts node ages to branch lengths

            Type    : Tree manipulator
            Title   : agetobl
            Usage   : $tree->agetobl;
            Function: Converts node ages to branch lengths
            Returns : The modified invocant.
            Args    : NONE
            Comments: This method uses ages as assigned to the generic 'age' slot
                      on the nodes in the trees. I.e. for each node in the tree,
                      $node->get_generic('age') must return a number

       ultrametricize()
           Sets all root-to-tip path lengths equal.

            Type    : Tree manipulator
            Title   : ultrametricize
            Usage   : $tree->ultrametricize;
            Function: Sets all root-to-tip path
                      lengths equal by stretching
                      all terminal branches to the
                      height of the tallest node.
            Returns : The modified invocant.
            Args    : NONE
            Comments: This method is analogous to
                      the 'ultrametricize' command
                      in Mesquite, i.e. no rate smoothing
                      or anything like that happens, just
                      a lengthening of terminal branches.

       scale()
           Scales the tree to the specified height.

            Type    : Tree manipulator
            Title   : scale
            Usage   : $tree->scale($height);
            Function: Scales the tree to the
                      specified height.
            Returns : The modified invocant.
            Args    : $height = a numerical value
                      indicating root-to-tip path length.
            Comments: This method uses the
                      $tree->calc_tree_height method, and
                      so for additive trees the *average*
                      root-to-tip path length is scaled to
                      $height (i.e. some nodes might be
                      taller than $height, others shorter).

       resolve()
           Randomly breaks polytomies.

            Type    : Tree manipulator
            Title   : resolve
            Usage   : $tree->resolve;
            Function: Randomly breaks polytomies by inserting
                      additional internal nodes.
            Returns : The modified invocant.
            Args    :
            Comments:

       prune_tips()
           Prunes argument nodes from invocant.

            Type    : Tree manipulator
            Title   : prune_tips
            Usage   : $tree->prune_tips(\@taxa);
            Function: Prunes specified taxa from invocant.
            Returns : A pruned Bio::Phylo::Forest::Tree object.
            Args    : A reference to an array of taxon names, or a taxa block, or a
                      reference to an array of taxon objects, or a reference to an
                      array of node objects
            Comments:

       keep_tips()
           Keeps argument nodes from invocant (i.e. prunes all others).

            Type    : Tree manipulator
            Title   : keep_tips
            Usage   : $tree->keep_tips(\@taxa);
            Function: Keeps specified taxa from invocant.
            Returns : The pruned Bio::Phylo::Forest::Tree object.
            Args    : Same as prune_tips, but with inverted meaning
            Comments:

       negative_to_zero()
           Converts negative branch lengths to zero.

            Type    : Tree manipulator
            Title   : negative_to_zero
            Usage   : $tree->negative_to_zero;
            Function: Converts negative branch
                      lengths to zero.
            Returns : The modified invocant.
            Args    : NONE
            Comments:

       ladderize()
           Sorts nodes in ascending (or descending) order of number of
           children.

            Type    : Tree manipulator
            Title   : ladderize
            Usage   : $tree->ladderize(1);
            Function: Sorts nodes
            Returns : The modified invocant.
            Args    : Optional, a true value to reverse the sort order

       sort_tips()
           Sorts nodes in (an approximation of) the provided ordering. Given
           an array reference of taxa, an array reference of name strings or a
           taxa object, this method attempts to order the tips in the same
           way. It does this by recursively computing the rank for all
           internal nodes by taking the average rank of its children. This
           results in the following orderings:

            (a,b,c,d,e,f); => $tree->sort_tips( [ qw(a c b f d e) ] ) => (a,c,b,f,d,e);

            (a,b,(c,d),e,f); => $tree->sort_tips( [ qw(a b e d c f) ] ); => (a,b,(e,(d,c)),f);

            ((a,b),((c,d),e),f); => $tree->sort_tips( [ qw(a e d c b f) ] ); => ((e,(d,c)),(a,b),f);

            Type    : Tree manipulator
            Title   : sort_tips
            Usage   : $tree->sort_tips($ordering);
            Function: Sorts nodes
            Returns : The modified invocant.
            Args    : Required, an array reference (or taxa object) whose ordering to match

       exponentiate()
           Raises branch lengths to argument.

            Type    : Tree manipulator
            Title   : exponentiate
            Usage   : $tree->exponentiate($power);
            Function: Raises branch lengths to $power.
            Returns : The modified invocant.
            Args    : A $power in any of perl's number formats.

       log_transform()
           Log argument base transform branch lengths.

            Type    : Tree manipulator
            Title   : log_transform
            Usage   : $tree->log_transform($base);
            Function: Log $base transforms branch lengths.
            Returns : The modified invocant.
            Args    : A $base in any of perl's number formats.

       remove_unbranched_internals()
           Collapses internal nodes with fewer than 2 children.

            Type    : Tree manipulator
            Title   : remove_unbranched_internals
            Usage   : $tree->remove_unbranched_internals;
            Function: Collapses internal nodes
                      with fewer than 2 children.
            Returns : The modified invocant.
            Args    : NONE
            Comments:

       remove_orphans()
           Removes all unconnected nodes.

            Type    : Tree manipulator
            Title   : remove_orphans
            Usage   : $tree->remove_orphans;
            Function: Removes all unconnected nodes
            Returns : The modified invocant.
            Args    : NONE
            Comments:

       deroot()
           Collapses one of the children of a basal bifurcation

            Type    : Tree manipulator
            Title   : deroot
            Usage   : $tree->deroot;
            Function: Removes root
            Returns : The modified invocant.
            Args    : Optional: node to collapse
            Comments:

   UTILITY METHODS
       clone()
           Clones invocant.

            Type    : Utility method
            Title   : clone
            Usage   : my $clone = $object->clone;
            Function: Creates a copy of the invocant object.
            Returns : A copy of the invocant.
            Args    : Optional: a hash of code references to
                      override reflection-based getter/setter copying

                      my $clone = $object->clone(
                          'set_forest' => sub {
                              my ( $self, $clone ) = @_;
                              for my $forest ( @{ $self->get_forests } ) {
                                  $clone->set_forest( $forest );
                              }
                          },
                          'set_matrix' => sub {
                              my ( $self, $clone ) = @_;
                              for my $matrix ( @{ $self->get_matrices } ) {
                                  $clone->set_matrix( $matrix );
                              }
                      );

            Comments: Cloning is currently experimental, use with caution.
                      It works on the assumption that the output of get_foo
                      called on the invocant is to be provided as argument
                      to set_foo on the clone - such as
                      $clone->set_name( $self->get_name ). Sometimes this
                      doesn't work, for example where this symmetry doesn't
                      exist, or where the return value of get_foo isn't valid
                      input for set_foo. If such a copy fails, a warning is
                      emitted. To make sure all relevant attributes are copied
                      into the clone, additional code references can be
                      provided, as in the example above. Typically, this is
                      done by overrides of this method in child classes.

   SERIALIZERS
       to_nexus()
           Serializes invocant to nexus string.

            Type    : Stringifier
            Title   : to_nexus
            Usage   : my $string = $tree->to_nexus;
            Function: Turns the invocant tree object
                      into a nexus string
            Returns : SCALAR
            Args    : Any arguments that can be passed to Bio::Phylo::Forest::to_nexus

       to_newick()
           Serializes invocant to newick string.

            Type    : Stringifier
            Title   : to_newick
            Usage   : my $string = $tree->to_newick;
            Function: Turns the invocant tree object
                      into a newick string
            Returns : SCALAR
            Args    : NONE

       to_xml()
           Serializes invocant to xml.

            Type    : Serializer
            Title   : to_xml
            Usage   : my $xml = $obj->to_xml;
            Function: Turns the invocant object into an XML string.
            Returns : SCALAR
            Args    : NONE

       to_svg()
           Serializes invocant to SVG.

            Type    : Serializer
            Title   : to_svg
            Usage   : my $svg = $obj->to_svg;
            Function: Turns the invocant object into an SVG string.
            Returns : SCALAR
            Args    : Same args as the Bio::Phylo::Treedrawer constructor
            Notes   : This will only work if you have the SVG module
                      from CPAN installed on your system.

       to_dom()
            Type    : Serializer
            Title   : to_dom
            Usage   : $tree->to_dom($dom)
            Function: Generates a DOM subtree from the invocant
                      and its contained objects
            Returns : an Element object
            Args    : DOM factory object

SEE ALSO
       There is a mailing list at
       https://groups.google.com/forum/#!forum/bio-phylo
       <https://groups.google.com/forum/#!forum/bio-phylo> for any user or
       developer questions and discussions.

       Bio::Phylo::Listable
           The Bio::Phylo::Forest::Tree object inherits from the
           Bio::Phylo::Listable object, so the methods defined therein also
           apply to trees.

       Bio::Phylo::Manual
           Also see the manual: Bio::Phylo::Manual and
           <http://rutgervos.blogspot.com>.

CITATION
       If you use Bio::Phylo in published research, please cite it:

       Rutger A Vos, Jason Caravas, Klaas Hartmann, Mark A Jensen and Chase
       Miller, 2011. Bio::Phylo - phyloinformatic analysis using Perl.  BMC
       Bioinformatics 12:63.  http://dx.doi.org/10.1186/1471-2105-12-63
       <http://dx.doi.org/10.1186/1471-2105-12-63>

perl v5.16.3                      2014-02-08   Bio::Phylo::Forest::TreeRole(3)
