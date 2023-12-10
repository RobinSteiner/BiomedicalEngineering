preamble
s = "The quick brown fox jumps over the lazy dog every Thursday and Saturday"
words = strsplit(s, ' ');
wordLengths = strlength(words);
shortestWords = words(wordLengths==min(wordLengths));
longestWords = words(wordLengths==max(wordLengths));
s1 = shortestWords(end)
s2 = longestWords(1)
s3 = join([upper(words(1)), words(2:end-1), upper(words(end))], ' ')
s4 = sort(words(lower(extract(words,1)) == extract(words,1)))
%%
preamble
N = 10000;
n = 4;
i = n:2:N;
prime_sums = primes(N) + primes(N)';
unique_sums = unique(prime_sums(mod(prime_sums, 2) == 0));
assert(all(i'==unique_sums((unique_sums>=n)&(unique_sums<=N))));