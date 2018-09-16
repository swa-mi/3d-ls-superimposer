# Superimposition code for cliques

Performs 3d least square superimposition of two sets of points of equal size (in terms of number of points) by translation and rotation, in constant time. Uses Simon K Kearsley's algorithm<sup>1</sup> for the purpose.

The code in `AlignmentFunctions.cpp` and `Align_cliques.cpp` may be of use to someone who wants a general implementation of Kearsley's algorithm. The rest of the code is for my own usage.

    Allowed options:
    -h [ --help ]               produce help message
    -i [ --input-cliqfile ] arg path to .cliq file
    -q [ --query-dir ] arg      Path to dir where query should be extracted from
    -t [ --target-dir ] arg     Path to dir where targets should be extracted 
    from
    -n [ --size ] arg (=8)      Set size of clique
    -r [ --rmsd ] arg (=2)      Set rmsd cutoff for succesful 
    superimposition/match
    -p [ --penalty ] arg (=10)  Set penalty for non-superimposition

Requires boost libraries.

Reads query from `--query-dir` argument, and superimposes that query set of points against targets in `--target-dir` argument
To test, run`cat test-files/20000_2mta_r3_r1_r1_r1_r2_r2_r2_r8.cliq.cliq test-files/20_2mta_r3_r1_r1_r1_r2_r2_r2_r8.cliq| bin/superimposer.bin -q test-files -t test-files` and compare the output against `test-files/expected_on20000.clique`

References:
1. On the orthogonal transformation used for structural comparisons. Kearsley, S. K., Acta Crystallographica Section A , Volume 45 (2) – Feb 1, 1989. doi:10.1107/s0108767388010128

Swastik Mishra, Parichit Sharma, Neelesh Soni
