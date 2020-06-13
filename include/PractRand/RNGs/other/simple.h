
/*
RNGs in the mediocre directory are not intended for real world use
only for research; as such they may get pretty sloppy in some areas

This set is of RNGs that do not make any significant use of:
	multiplication/division, arrays, flow control, complex math functions
*/

namespace PractRand {
	namespace RNGs {
		namespace Polymorphic {
			namespace NotRecommended {
				class xsalta16x3 : public vRNG16 {
					Uint16 a, b, c;
				public:
					Uint16 raw16();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class xsaltb16x3 : public vRNG16 {
					Uint16 a, b, c;
				public:
					Uint16 raw16();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class xsaltc16x3 : public vRNG16 {
					Uint16 a, b, c;
				public:
					Uint16 raw16();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				//xorshift RNGs, a subset of LFSRs proposed by Marsaglia in 2003
				class xorshift32 : public vRNG32 {
					//constants those Marsaglia described as "one of my favorites" on page 4 of his 2003 paper
					Uint32 a;
				public:
					Uint32 raw32();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class xorshift64 : public vRNG64 {
					//constants those Marsaglia used in his sample code on page 4 of his 2003 paper
					Uint64 a;
				public:
					Uint64 raw64();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class xorshift64of128 : public vRNG64 {
					//the constants are still in need of tuning
					Uint64 high, low;
					void xls(int bits);
					void xrs(int bits);
				public:
					Uint64 raw64();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class xorshift32of128 : public vRNG32 {
					xorshift64of128 impl;
				public:
					Uint32 raw32() {return Uint32(impl.raw64());}
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class xorshift16of32 : public vRNG16 {
					xorshift32 impl;
				public:
					Uint16 raw16() {return Uint16(impl.raw32());}
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class xorshift32of64 : public vRNG32 {
					xorshift64 impl;
				public:
					Uint32 raw32() {return Uint32(impl.raw64());}
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};

				class xorshift32x4 : public vRNG32 {
					//recommended at the top of Marsaglias 2003 xorshift paper
					Uint32 x,y,z,w;
				public:
					Uint32 raw32();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class xorwow32of96 : public vRNG32 {
					xorshift64 impl;
					Uint32 a;
				public:
					Uint32 raw32();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class xorwow32x6 : public vRNG32 {
					//recommended at the top of Marsaglias 2003 xorshift paper
					Uint32 x,y,z,w,v,d;
				public:
					Uint32 raw32();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};


				//by Ilya O. Levin, see http://www.literatecode.com/2004/10/18/sapparot/
				class sapparot : public vRNG32 {
					Uint32 a, b;
				public:
					Uint32 raw32();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};

				//variants of sapparot created for testing purposes
				class sap16of48 : public vRNG16 {
					Uint16 a, b, c;
				public:
					Uint16 raw16();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class sap32of96 : public vRNG32 {
					Uint32 a, b, c;
				public:
					Uint32 raw32();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				
				//the low quality variant of the FLEA RNG by Robert Jenkins
				//(he named this variant flea0)
				class flea32x1 : public vRNG32 {
					enum {SIZE = 1};
					Uint32 a[SIZE], b, c, d, i;
				public:
					Uint32 raw32();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};

				//versions of my SFC algorithm, archived here
				class sfc_v1_16 : public vRNG16 {
					Uint16 a, b, counter;
				public:
					Uint16 raw16();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class sfc_v1_32: public vRNG32 {
					Uint32 a, b, counter;
				public:
					Uint32 raw32();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class sfc_v2_16 : public vRNG16 {
					Uint16 a, b, counter;
				public:
					Uint16 raw16();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class sfc_v2_32 : public vRNG32 {
					Uint32 a, b, counter;
				public:
					Uint32 raw32();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class sfc_v3_16 : public vRNG16 {
					Uint16 a, b, counter;
				public:
					Uint16 raw16();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class sfc_v3_32 : public vRNG32 {
					Uint32 a, b, counter;
				public:
					Uint32 raw32();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class jsf16 : public vRNG16 {
					Uint16 a, b, c, d;
				public:
					Uint16 raw16();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};

				//a few simple RNGs just for testing purposes
				class simpleA : public vRNG32 {
					Uint32 a, b, c;
				public:
					Uint32 raw32();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class simpleB : public vRNG16 {
					Uint16 a, b, c;
				public:
					Uint16 raw16();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class simpleC : public vRNG16 {
					Uint16 a, b, c;
				public:
					Uint16 raw16();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class simpleD : public vRNG32 {
					Uint32 a, b, c;
				public:
					Uint32 raw32();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class simpleE : public vRNG32 {
					//seems like a good combination of speed & quality
					//but falls flat when used on 16 bit words (irreversible, statespace issues)
					Uint32 a, b, c;
				public:
					Uint32 raw32();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class simpleF : public vRNG16 {
					Uint16 a, b, c, d;
				public:
					Uint16 raw16();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class simpleG : public vRNG32 {
					Uint32 a, b, c, d;
				public:
					Uint32 raw32();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class trivium_weakenedA : public vRNG32 {
					Uint64 a, b;
				public:
					Uint32 raw32();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class trivium_weakenedB : public vRNG16 {
					Uint32 a, b, c;
				public:
					Uint16 raw16();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				//see http://www.drdobbs.com/tools/229625477
				class mo_Lesr32 : public vRNG32 {
					Uint32 state;
				public:
					Uint32 raw32();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class mo_ResrRers32 : public vRNG32 {
					Uint32 a;
					Uint32 b;
				public:
					Uint32 raw32();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class mo_Rers32of64 : public vRNG32 {
					Uint64 state;
				public:
					Uint32 raw32();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class mo_Resr32of64 : public vRNG32 {
					Uint64 state;
				public:
					Uint32 raw32();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class mo_Resdra32of64 : public vRNG32 {
					Uint64 state;
				public:
					Uint32 raw32();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
			}
		}
	}
}
