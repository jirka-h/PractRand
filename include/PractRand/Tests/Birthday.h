namespace PractRand {
	namespace Tests {

/*
Birthday tests
	originally I had little respect for the Birthday Test, as it looks kind of awful, and most implementations (on 32 bit integers) are kind of awful
		however, on finding a paper (I think it was by the TestU01 authors) suggesting that it works best at lambda~=1, I became more interested
			the results in that paper looked kind of awesome, and testing bore that out somewhat when using large sort buffers
			unfortunately, large sort buffers are... problematic
				PractRand perfers to assume cycles (and to a lesser extent, cache footprint) are the constraining factor in testing
					when users don't have fast enough CPUs, they can simply run tests longer ; when quality isn't high enough, they can simply run tests MUCH longer
				when the constraining factor is memory, scaling tests becomes impossible
				also, PractRand makes some choices that don't go well with memory constraints
					for instance, the "foldings" means that PractRand typically has 3 to 4 copies of each test in memory simultaneously
	lambda=1 (expected number of duplicate deltas of about 1.0... or at least keeping it between 0.5 and 4.0) means roughly a buffer size proportional to the cube root of the maximum word value
	something like a 2^11 word buffer for 32 bit words, 2^21 word buffer for 64 bit words, 2^43 word buffer for 128 bit words, etc
		which is obviously impractical for 128 bit words - that would be 128 terabytes assuming an in-place sorting algorithm was used
	intermediate word sizes are fine - just mask out some bits ; no need to fix up the alignments or anything, ignoring some bits is fine
	quality of test seems to scale quite decently with sort-buffer size, so long as word sizes are adjusted to keep lambda at least faintly near 1
	the actual test involves sorting the buffers, looking at the differences between adjacent sorted elements, 
		and then sorting the differences to check for duplicates
	BirthdaySystematic adds the minor innovation of permitting preliminary results to accessed even before the buffer is first filled
		this is done by masking out extra bits, so the word size is effectively smaller
		additionally, this means that the buffer is partially sorted afterwards, so ideally the later sortings should be faster
		and this means that the buffer has plenty of unused space in it, which the deltas can be stored in temporarily
		provided that preliminary results aren't needed after the buffer is more than half full, but before it's completed
	additional innovations are required
		since masking out bits and then looking for duplicates seems wasteful, I thought I'd leave as many bits in as possible
			and looks at the spacings between deltas instead of the duplicates of deltas (log and exp of delta between deltas)
			that's what I did in BirthdayAlt, but unfortunately, this did not seem to work well at all
			could try it by looking at the distrubtion of the number of bits that match between two adjacent sorted deltas
			but I was disappointed enough with BirthdayAlternative that I haven't felt up to it
		look for duplicates between different sets of deltas
			not entirely clear that this would be helpful, even if it does kinda approximate a larger sort buffer
			and the sizes of the sets of deltas are too large for this to be any more practical that larger sort buffers without alternative means of storing the deltas
			I considered storing only deltas that had shown up as duplicates, and looking for more matches of those
				but in preliminary testing it didn't look very good
			average delta value should be (2^word_bits / buffer_size)
				that puts it at something like 2^21 for 32 bit words, 2^43 for 64 bit words, and 2^85 for 128 bit words
				if the distribution is fairly tight then, instead of making a list of deltas and sorting them to look for duplicates
					we could make a bitmask of them (maybe plus a list of outliers)
						though it would have to be a VERY narrow distribution at 2^85 or even 2^43
						probably not practical for useful word sizes under normal circumstances
					we could make a hashtable of them... actually... that sounds promising
						index table by, say, bits 0-15 of deltas, fill table with bits 16-79 of deltas
							assuming expected deltas are large, and no weird structures in the fine details...
							it would probably be worse than simply using a large sort buffer
							but it might scale better when looking for matching deltas between multiple sort buffers
							or 
						number of times when overwriting the old value with the new changes nothing is the result
		filter input - only buffer lowest 0.1% of values (those that would end up at the start of the sorted buffer, so we can actually calculate the deltas we would have seen correctly... some of them anyway)
			possibly even increase the filtering over time, after we reach lambda=1 a few times or something
				that would let us mask out fewer bits over time while keeping the same buffer size
			would it actually help anything though?
				we can adjust bits discarded to keep lambda=1
				adjust filters to keep the buffer size constant as fewer bits are discarded on the low end
				but fundamentally, is quality tied to the buffer size, the number of deltas we can compare?
			possibly run at multiple filtering levels... actually that does sound good, provided that they're not too close
				say, one at 2^8, one at 2^16, one at 2^24, and a final one at 2^32
					going as high as 2^32 means that there would be useful new results coming in no matter how far
				an unfiltered version may or may not be useful
					skipping that would improve speed a lot, but unfiltered is the most useful at least early on
					and we're not actually sure filtered versions are any good yet
				could handle the higher filter versions only once the lower filter versions are sorted
					that has the advantage of being faster, and minor disadvantage of delaying higher level results
					... but it's probably kind of pointless, given that higher filterings are so rarely hit anyway
					would probably only be worthwhile if the filtering levels weren't very different
					but in that case, we'd run out of memory from having too many sort buffers
*/
		class Birthday32 : public TestBaseclass {
			enum {
				BUFFER_SIZE_L2 = 12, // must be at least 8
				BUFFER_SIZE = 1 << BUFFER_SIZE_L2,
				MAX_DUPLICATES = 32
			};
			Uint32 buffer[1 << BUFFER_SIZE_L2];
			Uint64 counts[MAX_DUPLICATES];
			int num_buffered;
			void flush_buffer();
		public:
			Birthday32();
			virtual void init(PractRand::RNGs::vRNG *known_good);
			//virtual void deinit();
			virtual std::string get_name() const;
			virtual void get_results(std::vector<TestResult> &results);

			virtual void test_blocks(TestBlock *data, int numblocks);
		};
		class Birthday64 : public TestBaseclass {
			enum {
				BUFFER_SIZE_L2 = 23, // must be at least 7
				BUFFER_SIZE = 1 << BUFFER_SIZE_L2,
				MAX_DUPLICATES = 64,
				SORT_HELPER_BITS = 10
			};
			Uint64 buffer[1 << BUFFER_SIZE_L2];
			Uint64 counts[MAX_DUPLICATES];
			static void _histogram_in_place_sort64(Uint64 *base, long length, long bits_already, Uint32 freq_counts[1 << SORT_HELPER_BITS]);
			static void _histogram_in_place_sort64(Uint64 *base, long length);
			static void _histogram_sort64(Uint64 *base, long length, long bits_already, Uint32 freq_counts[1 << SORT_HELPER_BITS]);
			static void _histogram_sort64(Uint64 *base, long length);
			int num_buffered;
			void flush_buffer();
		public:
			Birthday64();
			virtual void init(PractRand::RNGs::vRNG *known_good);
			//virtual void deinit();
			virtual std::string get_name() const;
			virtual void get_results(std::vector<TestResult> &results);

			virtual void test_blocks(TestBlock *data, int numblocks);
		};
		namespace BirthdayHelpers {
			enum { SORT_HELPER_BITS = 8 };
			struct i128 {
				Uint64 low;
				Uint64 high;
				bool operator==(const i128 &other) const {
					return high == other.high && low == other.low;
				}
				bool operator<(const i128 &other) const {
					if (high < other.high) return true;
					if (high > other.high) return false;
					return low < other.low;
				}
				i128 operator-(const i128 &other) const {
					i128 rv;
					rv.high = high - other.high;
					rv.low = low - other.low;
					rv.high -= other.low > low ? 1 : 0;
					return rv;
				}
			};

			// this is the fastest in-place sort I've tried so far
			// in-place is important if I want to sort something huge without assuming I can allocate a comparable amount of memory to help with the sorting
			void histogram_in_place_sort128(i128 *base, Uint64 length, long bits_already, Uint64 freq_counts[1 << SORT_HELPER_BITS]);
			void histogram_in_place_sort128(i128 *base, Uint64 length, long bits_already = 0);

			// ..but sometimes I'm sorting smaller buffers, with pre-allocated regions to sort into..., so maybe another interface for that would help
			void histogram_sort_and_copy(i128 *base, i128 *dest, Uint64 length, long bits_already, Uint64 freq_counts[1 << SORT_HELPER_BITS]);
			void histogram_sort_and_copy(i128 *base, i128 *dest, Uint64 length, long bits_already = 0);
			//possibly faster algorithm for the same interface?
			void radix_sort_and_copy(i128 *base, i128 *dest, Uint64 length, long bits_already = 0);

			void _sorted_deltas_of_sorted_values(i128 *base, long length_L2, Uint64 freq_counts[1 << SORT_HELPER_BITS]);
			void _sorted_deltas_of_sorted_values(i128 *base, long length_L2);
		};
		class BirthdayLamda1 : public TestBaseclass {
		protected:
			//optimized for lambda=1, few runs, as described in "On the performance of birthday spacings tests with certain families of random number generators" (L'ecuyer & Simard, 2001)
			typedef BirthdayHelpers::i128 i128;
			enum {
				SORT_HELPER_BITS = BirthdayHelpers::SORT_HELPER_BITS,
				DO_LARGEST_SPACING = 1,
			};
			bool autofail;
			Uint64 sort_helper_counts[1 << SORT_HELPER_BITS];
			//i128 buffer[1 << BUFFER_SIZE_L2];//can't have arrays this large inside a class due to object file or executable file format constraints
			//std::vector<i128> buffer;// ... and the STL vector implementation I'm using throws some kind of exception if it exceeds about 4 GB or so
			i128 *buffer;
			Uint64 num_buffered;
			virtual Uint64 flush_buffer();
			double duplicates;
			double expected_duplicates;
			double longest_spacing;
			int buffer_size_L2;
			int bits_to_use;
		public:
			BirthdayLamda1(int buffer_size_L2_ = 26);
			virtual ~BirthdayLamda1();
			virtual void init(PractRand::RNGs::vRNG *known_good);
			//virtual void deinit();
			virtual std::string get_name() const;
			virtual void get_results(std::vector<TestResult> &results);

			virtual void test_blocks(TestBlock *data, int numblocks);
		};
		/*
		class BirthdaySystematic64 {
			// does Birthday Spacings test on 64 bit integers
			// if a result is requested mid-buffer, the buffer will be sorted and processed but the contents left reusable and the results treated as ephemeral
			// generally, the buffer size is aimed to have a lambda value of about 1 when using all 64 bits
			// multiple tests are combined by simply adding their observed duplicated and expected duplicates
			// currently undecided on whether or not early use of the buffer will suppress some bits or not
			enum {
				BUFSIZE_L2 = 22, 
				BUFSIZE = 1 << BUFSIZE_L2
			};
			Uint64 buffer[BUFSIZE];
			Uint64 elements_buffered;
			Uint64 num_sorted; // this many elements in the buffer are already sorted, starting at the beginning, potentially allowing optimization to the final sorting of its contents
			Uint64 observed_duplicates;
			double expected_duplicates;
			Uint64 evaluate_buffer();
		public:
			;
		};*/
		class BirthdaySystematic128 : public BirthdayLamda1 {
			// similar to BirthdayLambda1 above
			// but if a result is requested before the first sample is ready, it will return a result for a partial buffer
			// and attempts to have everything optimized for the possibility of that partial-buffer case
			virtual Uint64 flush_buffer();
			static Uint64 get_target_num_at_bufsize(int bufsize_L2_);
			int already_sorted;//if this is half of (1ull << bufsize_L2) then incomplete_duplicates should hold 

			double score;//for scoring method 2
			static double evaluate_score(double lambda, Uint64 duplicates);

			void do_incomplete_buffer();
			double incomplete_duplicates;
			double incomplete_expected_duplicates;
		public:
			BirthdaySystematic128(int max_bufsize_L2_ = 28);
			virtual void init(PractRand::RNGs::vRNG *known_good);
			virtual std::string get_name() const;
			virtual void get_results(std::vector<TestResult> &results);
			virtual void test_blocks(TestBlock *data, int numblocks);
		};
		class BirthdayAlt : public TestBaseclass {
			//as for BirthdayLambda1, but: 
			// keep all bits regardless of buffer size, just count an adjusting range of near deltas as if they were exact matches (or score them based upon how exact they are?)
			// try filtering the initial samples range, as if it was a small part of a larger sort buffer
			typedef BirthdayHelpers::i128 i128;
			enum {
				SORT_HELPER_BITS = BirthdayHelpers::SORT_HELPER_BITS,
				//FILTER_BITS = 0,
			};
			//i128 buffer[1 << BUFFER_SIZE_L2];
			//std::vector<i128> buffer;
			i128 *buffer;
			int num_buffered;
			int buffer_size_L2;
			int filter_bits;
			Uint64 sort_helper_counts[1 << SORT_HELPER_BITS];
			bool autofail;
			void flush_buffer();

			double score_sum_log;
			double score_sum_log2;
			double score_sum_log_sqr;
			Uint64 count;
			static void _lookup_constants(int buffer_size_L2, long double *offset, long double *deviation, long double *samples);
		public:
			BirthdayAlt(int buffer_size_L2_, int filter_bits_ = 0);
			~BirthdayAlt();
			virtual void init(PractRand::RNGs::vRNG *known_good);
			//virtual void deinit();
			virtual std::string get_name() const;
			virtual void get_results(std::vector<TestResult> &results);

			virtual void test_blocks(TestBlock *data, int numblocks);
		};//*/
	}//Tests
}//PractRand
