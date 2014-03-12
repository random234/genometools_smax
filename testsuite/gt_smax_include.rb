Name "gt smax map only index"
Keywords "gt_smax"
Test do 
  run_test "#{$bin}gt suffixerator -db #{$testdata}Random.fna -dna -suf -sds" +
    " -lcp -tis -des -ssp"   
  run_test "#{$bin}gt smax -ii Random.fna -map"
end

Name "gt smax map absolute"
Keywords "gt_smax"
Test do 
  run_test "#{$bin}gt suffixerator -db #{$testdata}Random.fna -dna -suf -sds" +
    " -lcp -tis -des -ssp"   
  run_test "#{$bin}gt smax -ii Random.fna -absolute -map"
end

Name "gt smax map diff absolute"
Keywords "gt_smax"
Test do
  run_test "#{$bin}gt suffixerator -db #{$testdata}Random.fna -dna -suf -sds" +
    " -lcp -tis -des -ssp"
  run_test "#{$bin}gt smax -ii Random.fna -absolute -l 2 -map"
  run_test "diff #{last_stdout} #{$testdata}Random.fna_smax_absolute.out"
end

Name "gt smax map diff relative"
Keywords "gt_smax"
Test do
  run_test "#{$bin}gt suffixerator -db #{$testdata}Random.fna -dna -suf -sds" +
    " -lcp -tis -des -ssp"
  run_test "#{$bin}gt smax -ii Random.fna -l 2 -map"
  run_test "diff #{last_stdout} #{$testdata}Random.fna_smax_relative.out"
end
