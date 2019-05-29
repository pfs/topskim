#ifndef tnp_weight_h
#define tnp_weight_h

#include "TMath.h"

// IN THIS FILE YOU WILL FIND:
// ++++++++++++++
//
// - MuID: (tnp_weight_muid_pbpb)   Preliminary
//   * idx = 0: nominal
//   * idx = -1: syst variation,  +1 sigma
//   * idx = -2: syst variation,  -1 sigma
//   * idx = +1: stat variation,  +1 sigma
//   * idx = +2: stat variation,  -1 sigma
//
// - Trigger: (tnp_weight_trg_pbpb)   NOT UPDATED
//   * idx = 0:  nominal
//   * idx = 1..100: toy variations, stat. only
//   * idx = -1: syst variation, "new_MAX", +1 sigma
//   * idx = -2: syst variation, "new_MAX", -1 sigma
//   * idx = -10: binned

// For all:
//   * idx = +200: tnp efficiency for data
//   * idx = +300: tnp efficiency for MC

// THE INDIVIDUAL SFs
// ++++++++++++++++++
double tnp_weight_muid_pbpb(double eta, int idx=0);
double tnp_weight_trg_pbpb(double pt, double eta, int idx=0);  // NOT UPDATED



///////////////////////////////////////////////////
//                 M u I D    P b P b                //
///////////////////////////////////////////////////
double tnp_weight_muid_pbpb(double eta, int idx)
{
   double x = eta;
   double syst = 0.6e-2;  //preliminary

   double num=1,den=1;
   

   // MC
   if (x > -2.4&&x <= -2.1) den = 0.994717;
   if (x > -2.1&&x <= -1.6) den = 0.993653;
   if (x > -1.6&&x <= -1.2) den = 0.982687;
   if (x > -1.2&&x <= -0.9) den = 0.962992;
   if (x > -0.9&&x <= -0.6) den = 0.970103;
   if (x > -0.6&&x <= -0.3) den = 0.982073;
   if (x > -0.3&&x <= 0) den = 0.968219;
   if (x > 0 && x <= 0.3) den = 0.961703;
   if (x > 0.3&&x <= 0.6) den = 0.978947;
   if (x > 0.6&&x <= 0.9) den = 0.970324;
   if (x > 0.9&&x <= 1.2) den = 0.956433;
   if (x > 1.2&&x <= 1.6) den = 0.982577;
   if (x > 1.6&&x <= 2.1) den = 0.995146;
   if (x > 2.1&&x <= 2.4) den = 0.994106;


   // data
   if (idx <= 0 || idx > 10) { // nominal
	   if (x > -2.4&&x <= -2.1) num = 0.983734;
	   if (x > -2.1&&x <= -1.6) num = 0.993892;
	   if (x > -1.6&&x <= -1.2) num = 0.979222;
	   if (x > -1.2&&x <= -0.9) num = 0.955529;
	   if (x > -0.9&&x <= -0.6) num = 0.966578;
	   if (x > -0.6&&x <= -0.3) num = 0.983958;
	   if (x > -0.3&&x <= 0) num = 0.959491;
	   if (x > 0 && x <= 0.3) num = 0.959744;
	   if (x > 0.3&&x <= 0.6) num = 0.976402;
	   if (x > 0.6&&x <= 0.9) num = 0.96839;
	   if (x > 0.9&&x <= 1.2) num = 0.961177;
	   if (x > 1.2&&x <= 1.6) num = 0.98078;
	   if (x > 1.6&&x <= 2.1) num = 0.991852;
	   if (x > 2.1&&x <= 2.4) num = 0.993082;
   }
   else if (idx == 1) { // stat up
	   if (x > -2.4&&x <= -2.1) num = 0.98710825;
	   if (x > -2.1&&x <= -1.6) num = 0.99575701;
	   if (x > -1.6&&x <= -1.2) num = 0.98250043;
	   if (x > -1.2&&x <= -0.9) num = 0.96029298;
	   if (x > -0.9&&x <= -0.6) num = 0.97057735;
	   if (x > -0.6&&x <= -0.3) num = 0.98705592;
	   if (x > -0.3&&x <= 0) num = 0.96370511;
	   if (x > 0 && x <= 0.3) num = 0.96396652;
	   if (x > 0.3&&x <= 0.6) num = 0.97980996;
	   if (x > 0.6&&x <= 0.9) num = 0.97234624;
	   if (x > 0.9&&x <= 1.2) num = 0.9659164;
	   if (x > 1.2&&x <= 1.6) num = 0.98388187;
	   if (x > 1.6&&x <= 2.1) num = 0.99371692;
	   if (x > 2.1&&x <= 2.4) num = 0.9961665;
   }
   else if (idx == 2) { // stat down
	   if (x > -2.4&&x <= -2.1) num = 0.98075691;
	   if (x > -2.1&&x <= -1.6) num = 0.99224693;
	   if (x > -1.6&&x <= -1.2) num = 0.97618383;
	   if (x > -1.2&&x <= -0.9) num = 0.95101985;
	   if (x > -0.9&&x <= -0.6) num = 0.96279449;
	   if (x > -0.6&&x <= -0.3) num = 0.98108946;
	   if (x > -0.3&&x <= 0) num = 0.95550071;
	   if (x > 0 && x <= 0.3) num = 0.95575316;
	   if (x > 0.3&&x <= 0.6) num = 0.97322907;
	   if (x > 0.6&&x <= 0.9) num = 0.96468663;
	   if (x > 0.9&&x <= 1.2) num = 0.9567025;
	   if (x > 1.2&&x <= 1.6) num = 0.9779218;
	   if (x > 1.6&&x <= 2.1) num = 0.99018962;
	   if (x > 2.1&&x <= 2.4) num = 0.9905424;
   }

   if (idx == 200) den = 1.;
   if (idx == 300) num = den * den;


   double syst_factor = 1.;
   if (idx == -1) syst_factor = 1 + syst;
   if (idx == -2) syst_factor = 1 - syst;
   return (num / den)*syst_factor;
}


///////////////////////////////////////////////////
//               T R G      P b P b              //
///////////////////////////////////////////////////
double tnp_weight_trg_pbpb(double pt, double eta, int idx) {
   double x = pt;

   if (idx == -10) {
      if (fabs(eta)<1.2) {
         // 0 < |eta| < 1.2
         if (pt<15) return 0.981114;
         else if (pt<25) return 1.00118;
         else if (pt<30) return 0.99319;
         else if (pt<35) return 0.994709;
         else if (pt<40) return 0.991174;
         else if (pt<45) return 0.993454;
         else if (pt<50) return 0.987039;
         else if (pt<80) return 0.992419;
         else if (pt<200) return 0.971874;
      } else if (fabs(eta)<2.1) {
         // 1.2 < |eta| < 2.1
         if (pt<15) return 1.0036;
         else if (pt<25) return 0.998054;
         else if (pt<30) return 0.995244;
         else if (pt<35) return 0.992706;
         else if (pt<40) return 0.999179;
         else if (pt<45) return 0.999337;
         else if (pt<50) return 0.995187;
         else if (pt<80) return 1.00077;
         else if (pt<200) return 1.00301;
      } else {
         // 2.1 < |eta| < 2.4
         if (pt<15) return 1.00886;
         else if (pt<25) return 1.00277;
         else if (pt<40) return 0.9927;
         else if (pt<80) return 0.985097;
         else if (pt<200) return 0.974808;
      }
   }

   double num = 1, den = 1;

   // nominal
   if (fabs(eta)<1.2) {
      num = 0.98711-0.00022*x;
      den = 0.98864-0.00002*x;
   } else if (fabs(eta)<2.1) {
      num = 0.99296 + 0.00001*x;
      den = 0.99633-0.00004*x;
   } else {
      num = 1.00072-0.00052*x;
      den = 0.99254-0.00017*x;
   }

   // 100 variations 
   if (fabs(eta)<1.2) {
      if (idx == 1   ) num = 0.98216-0.00012*x;
      else if (idx == 2   ) num = 0.98444-0.00015*x;
      else if (idx == 3   ) num = 0.99079-0.00033*x;
      else if (idx == 4   ) num = 0.97209+0.00010*x;
      else if (idx == 5   ) num = 0.99083-0.00029*x;
      else if (idx == 6   ) num = 0.99184-0.00029*x;
      else if (idx == 7   ) num = 0.98340-0.00016*x;
      else if (idx == 8   ) num = 0.98541-0.00026*x;
      else if (idx == 9   ) num = 0.97642+0.00005*x;
      else if (idx == 10  ) num = 0.99088-0.00026*x;
      else if (idx == 11  ) num = 0.97731-0.00003*x;
      else if (idx == 12  ) num = 0.98927-0.00032*x;
      else if (idx == 13  ) num = 0.99025-0.00031*x;
      else if (idx == 14  ) num = 0.99398-0.00040*x;
      else if (idx == 15  ) num = 0.98014-0.00017*x;
      else if (idx == 16  ) num = 0.99707-0.00046*x;
      else if (idx == 17  ) num = 0.97607-0.00002*x;
      else if (idx == 18  ) num = 0.99917-0.00049*x;
      else if (idx == 19  ) num = 0.98614-0.00018*x;
      else if (idx == 20  ) num = 0.99061-0.00024*x;
      else if (idx == 21  ) num = 0.98719-0.00022*x;
      else if (idx == 22  ) num = 0.98611-0.00021*x;
      else if (idx == 23  ) num = 0.98541-0.00014*x;
      else if (idx == 24  ) num = 0.98480-0.00015*x;
      else if (idx == 25  ) num = 0.97988-0.00009*x;
      else if (idx == 26  ) num = 0.97945-0.00006*x;
      else if (idx == 27  ) num = 0.99505-0.00038*x;
      else if (idx == 28  ) num = 0.99933-0.00049*x;
      else if (idx == 29  ) num = 0.98003-0.00000*x;
      else if (idx == 30  ) num = 0.97219+0.00009*x;
      else if (idx == 31  ) num = 0.98237-0.00010*x;
      else if (idx == 32  ) num = 0.98512-0.00015*x;
      else if (idx == 33  ) num = 0.98142-0.00007*x;
      else if (idx == 34  ) num = 0.99604-0.00039*x;
      else if (idx == 35  ) num = 0.98885-0.00031*x;
      else if (idx == 36  ) num = 0.99388-0.00034*x;
      else if (idx == 37  ) num = 0.99168-0.00033*x;
      else if (idx == 38  ) num = 0.97524+0.00004*x;
      else if (idx == 39  ) num = 0.98027-0.00011*x;
      else if (idx == 40  ) num = 0.99009-0.00022*x;
      else if (idx == 41  ) num = 0.98929-0.00026*x;
      else if (idx == 42  ) num = 0.98244-0.00015*x;
      else if (idx == 43  ) num = 0.98258-0.00011*x;
      else if (idx == 44  ) num = 0.98014-0.00002*x;
      else if (idx == 45  ) num = 0.98656-0.00026*x;
      else if (idx == 46  ) num = 0.98013+0.00004*x;
      else if (idx == 47  ) num = 0.98895-0.00027*x;
      else if (idx == 48  ) num = 0.98689-0.00028*x;
      else if (idx == 49  ) num = 0.98938-0.00026*x;
      else if (idx == 50  ) num = 0.97786-0.00007*x;
      else if (idx == 51  ) num = 0.97637+0.00003*x;
      else if (idx == 52  ) num = 0.99635-0.00050*x;
      else if (idx == 53  ) num = 0.98092-0.00011*x;
      else if (idx == 54  ) num = 0.98849-0.00032*x;
      else if (idx == 55  ) num = 0.98958-0.00022*x;
      else if (idx == 56  ) num = 0.99087-0.00023*x;
      else if (idx == 57  ) num = 0.99112-0.00035*x;
      else if (idx == 58  ) num = 0.99542-0.00035*x;
      else if (idx == 59  ) num = 0.99380-0.00035*x;
      else if (idx == 60  ) num = 0.98629-0.00019*x;
      else if (idx == 61  ) num = 0.97954-0.00005*x;
      else if (idx == 62  ) num = 0.99353-0.00037*x;
      else if (idx == 63  ) num = 0.99171-0.00033*x;
      else if (idx == 64  ) num = 0.99339-0.00038*x;
      else if (idx == 65  ) num = 0.98256-0.00003*x;
      else if (idx == 66  ) num = 0.98016-0.00020*x;
      else if (idx == 67  ) num = 0.99609-0.00037*x;
      else if (idx == 68  ) num = 0.99209-0.00030*x;
      else if (idx == 69  ) num = 1.00039-0.00055*x;
      else if (idx == 70  ) num = 0.98237-0.00018*x;
      else if (idx == 71  ) num = 0.99425-0.00040*x;
      else if (idx == 72  ) num = 0.98840-0.00017*x;
      else if (idx == 73  ) num = 0.99310-0.00048*x;
      else if (idx == 74  ) num = 0.98728-0.00028*x;
      else if (idx == 75  ) num = 0.99716-0.00041*x;
      else if (idx == 76  ) num = 0.99419-0.00036*x;
      else if (idx == 77  ) num = 0.98171-0.00012*x;
      else if (idx == 78  ) num = 0.98239-0.00006*x;
      else if (idx == 79  ) num = 0.98080-0.00012*x;
      else if (idx == 80  ) num = 0.99226-0.00032*x;
      else if (idx == 81  ) num = 0.98672-0.00018*x;
      else if (idx == 82  ) num = 0.99198-0.00048*x;
      else if (idx == 83  ) num = 0.98282-0.00018*x;
      else if (idx == 84  ) num = 0.98647-0.00020*x;
      else if (idx == 85  ) num = 0.99965-0.00056*x;
      else if (idx == 86  ) num = 1.00148-0.00056*x;
      else if (idx == 87  ) num = 0.97968-0.00012*x;
      else if (idx == 88  ) num = 0.98417-0.00014*x;
      else if (idx == 89  ) num = 0.98236-0.00010*x;
      else if (idx == 90  ) num = 0.98229-0.00015*x;
      else if (idx == 91  ) num = 0.98788-0.00021*x;
      else if (idx == 92  ) num = 0.98183-0.00011*x;
      else if (idx == 93  ) num = 0.98130-0.00005*x;
      else if (idx == 94  ) num = 0.97616-0.00006*x;
      else if (idx == 95  ) num = 0.98670-0.00024*x;
      else if (idx == 96  ) num = 0.98200-0.00012*x;
      else if (idx == 97  ) num = 0.98162-0.00010*x;
      else if (idx == 98  ) num = 0.97650+0.00003*x;
      else if (idx == 99  ) num = 0.99554-0.00045*x;
      else if (idx == 100 ) num = 0.97934-0.00010*x;
      // // ADD THE FUNCTIONS FOR SYSTS BELOW! +1 SIGMA (IDX==-1) AND -1 SIGMA (IDX==-2)
      else if (idx == -1  ) num = 0.99313-0.00032*x;
      else if (idx == -2  ) num = 0.98101-0.00013*x;
   } else if (fabs(eta)<2.1) {
      if (idx == 1   ) num = 0.99035+0.00003*x;
      else if (idx == 2   ) num = 0.99224+0.00003*x;
      else if (idx == 3   ) num = 0.99557-0.00009*x;
      else if (idx == 4   ) num = 0.99262-0.00003*x;
      else if (idx == 5   ) num = 0.98790+0.00006*x;
      else if (idx == 6   ) num = 0.99470-0.00010*x;
      else if (idx == 7   ) num = 0.99488-0.00004*x;
      else if (idx == 8   ) num = 0.99265+0.00000*x;
      else if (idx == 9   ) num = 0.99252-0.00001*x;
      else if (idx == 10  ) num = 0.99290+0.00003*x;
      else if (idx == 11  ) num = 0.99594-0.00013*x;
      else if (idx == 12  ) num = 0.99790-0.00019*x;
      else if (idx == 13  ) num = 0.99244+0.00000*x;
      else if (idx == 14  ) num = 0.99034+0.00006*x;
      else if (idx == 15  ) num = 0.99067+0.00003*x;
      else if (idx == 16  ) num = 0.99296-0.00000*x;
      else if (idx == 17  ) num = 0.99447-0.00012*x;
      else if (idx == 18  ) num = 0.99822-0.00015*x;
      else if (idx == 19  ) num = 0.99034+0.00004*x;
      else if (idx == 20  ) num = 0.98905+0.00009*x;
      else if (idx == 21  ) num = 0.99728-0.00013*x;
      else if (idx == 22  ) num = 0.98986+0.00006*x;
      else if (idx == 23  ) num = 0.99346-0.00001*x;
      else if (idx == 24  ) num = 0.99091+0.00005*x;
      else if (idx == 25  ) num = 0.99422-0.00007*x;
      else if (idx == 26  ) num = 0.99192+0.00004*x;
      else if (idx == 27  ) num = 0.99295+0.00000*x;
      else if (idx == 28  ) num = 0.99333+0.00003*x;
      else if (idx == 29  ) num = 0.99803-0.00013*x;
      else if (idx == 30  ) num = 0.99531-0.00004*x;
      else if (idx == 31  ) num = 0.99849-0.00018*x;
      else if (idx == 32  ) num = 1.00133-0.00023*x;
      else if (idx == 33  ) num = 0.99100+0.00004*x;
      else if (idx == 34  ) num = 0.99552-0.00012*x;
      else if (idx == 35  ) num = 0.99358-0.00003*x;
      else if (idx == 36  ) num = 0.99728-0.00016*x;
      else if (idx == 37  ) num = 0.99082+0.00004*x;
      else if (idx == 38  ) num = 0.98951+0.00009*x;
      else if (idx == 39  ) num = 0.99028+0.00005*x;
      else if (idx == 40  ) num = 0.99038+0.00006*x;
      else if (idx == 41  ) num = 0.99066+0.00004*x;
      else if (idx == 42  ) num = 0.99353-0.00007*x;
      else if (idx == 43  ) num = 0.99252-0.00002*x;
      else if (idx == 44  ) num = 0.99347-0.00003*x;
      else if (idx == 45  ) num = 0.99327+0.00001*x;
      else if (idx == 46  ) num = 1.00178-0.00033*x;
      else if (idx == 47  ) num = 0.99974-0.00021*x;
      else if (idx == 48  ) num = 0.99517-0.00011*x;
      else if (idx == 49  ) num = 0.99390-0.00002*x;
      else if (idx == 50  ) num = 0.99083+0.00003*x;
      else if (idx == 51  ) num = 0.99079+0.00004*x;
      else if (idx == 52  ) num = 0.99611-0.00004*x;
      else if (idx == 53  ) num = 0.99067+0.00006*x;
      else if (idx == 54  ) num = 0.99370-0.00011*x;
      else if (idx == 55  ) num = 0.99157+0.00004*x;
      else if (idx == 56  ) num = 0.99434-0.00006*x;
      else if (idx == 57  ) num = 0.99419-0.00006*x;
      else if (idx == 58  ) num = 0.99454-0.00002*x;
      else if (idx == 59  ) num = 0.99427-0.00008*x;
      else if (idx == 60  ) num = 0.99417-0.00002*x;
      else if (idx == 61  ) num = 0.99198+0.00003*x;
      else if (idx == 62  ) num = 0.99311-0.00007*x;
      else if (idx == 63  ) num = 0.99782-0.00020*x;
      else if (idx == 64  ) num = 0.99705-0.00017*x;
      else if (idx == 65  ) num = 0.99078+0.00001*x;
      else if (idx == 66  ) num = 0.99604-0.00009*x;
      else if (idx == 67  ) num = 0.98999+0.00007*x;
      else if (idx == 68  ) num = 0.99561-0.00009*x;
      else if (idx == 69  ) num = 1.00041-0.00026*x;
      else if (idx == 70  ) num = 0.99434-0.00004*x;
      else if (idx == 71  ) num = 0.99150+0.00000*x;
      else if (idx == 72  ) num = 0.99620-0.00015*x;
      else if (idx == 73  ) num = 0.99408-0.00001*x;
      else if (idx == 74  ) num = 0.99139+0.00004*x;
      else if (idx == 75  ) num = 0.99662-0.00013*x;
      else if (idx == 76  ) num = 0.99270+0.00001*x;
      else if (idx == 77  ) num = 0.99150-0.00001*x;
      else if (idx == 78  ) num = 0.99159-0.00001*x;
      else if (idx == 79  ) num = 0.99518-0.00006*x;
      else if (idx == 80  ) num = 0.99239+0.00003*x;
      else if (idx == 81  ) num = 0.99147+0.00004*x;
      else if (idx == 82  ) num = 0.99207+0.00002*x;
      else if (idx == 83  ) num = 0.99147+0.00001*x;
      else if (idx == 84  ) num = 0.99450-0.00002*x;
      else if (idx == 85  ) num = 0.99665-0.00010*x;
      else if (idx == 86  ) num = 0.99218+0.00003*x;
      else if (idx == 87  ) num = 0.99033+0.00005*x;
      else if (idx == 88  ) num = 0.99448-0.00005*x;
      else if (idx == 89  ) num = 0.98874+0.00005*x;
      else if (idx == 90  ) num = 0.99019+0.00005*x;
      else if (idx == 91  ) num = 0.99692-0.00016*x;
      else if (idx == 92  ) num = 0.99135+0.00002*x;
      else if (idx == 93  ) num = 0.99239-0.00001*x;
      else if (idx == 94  ) num = 0.99123+0.00005*x;
      else if (idx == 95  ) num = 0.98927+0.00008*x;
      else if (idx == 96  ) num = 0.99349-0.00001*x;
      else if (idx == 97  ) num = 0.99073+0.00005*x;
      else if (idx == 98  ) num = 1.00333-0.00029*x;
      else if (idx == 99  ) num = 0.99071+0.00006*x;
      else if (idx == 100 ) num = 0.99103+0.00002*x;
      // // ADD THE FUNCTIONS FOR SYSTS BELOW! +1 SIGMA (IDX==-1) AND -1 SIGMA (IDX==-2)
      else if (idx == -1  ) num = 0.99234+0.00003*x;
      else if (idx == -2  ) num = 0.99346-0.00001*x;
   } else {
      if (idx == 1   ) num = 0.98883-0.00012*x;
      else if (idx == 2   ) num = 0.99973-0.00052*x;
      else if (idx == 3   ) num = 1.00163-0.00073*x;
      else if (idx == 4   ) num = 1.00401-0.00077*x;
      else if (idx == 5   ) num = 1.00779-0.00083*x;
      else if (idx == 6   ) num = 0.98365+0.00006*x;
      else if (idx == 7   ) num = 0.99932-0.00061*x;
      else if (idx == 8   ) num = 1.01102-0.00087*x;
      else if (idx == 9   ) num = 0.98429+0.00004*x;
      else if (idx == 10  ) num = 1.00947-0.00100*x;
      else if (idx == 11  ) num = 1.01235-0.00110*x;
      else if (idx == 12  ) num = 1.00550-0.00077*x;
      else if (idx == 13  ) num = 0.98864-0.00009*x;
      else if (idx == 14  ) num = 1.00527-0.00073*x;
      else if (idx == 15  ) num = 0.99135-0.00022*x;
      else if (idx == 16  ) num = 0.99817-0.00033*x;
      else if (idx == 17  ) num = 1.00889-0.00088*x;
      else if (idx == 18  ) num = 0.98720-0.00007*x;
      else if (idx == 19  ) num = 1.00314-0.00068*x;
      else if (idx == 20  ) num = 0.98765-0.00013*x;
      else if (idx == 21  ) num = 1.01230-0.00097*x;
      else if (idx == 22  ) num = 1.01536-0.00112*x;
      else if (idx == 23  ) num = 0.98338+0.00002*x;
      else if (idx == 24  ) num = 0.99702-0.00047*x;
      else if (idx == 25  ) num = 1.00064-0.00067*x;
      else if (idx == 26  ) num = 0.99694-0.00036*x;
      else if (idx == 27  ) num = 1.01886-0.00139*x;
      else if (idx == 28  ) num = 1.01914-0.00120*x;
      else if (idx == 29  ) num = 0.99698-0.00036*x;
      else if (idx == 30  ) num = 0.99037-0.00022*x;
      else if (idx == 31  ) num = 0.99771-0.00039*x;
      else if (idx == 32  ) num = 1.00317-0.00080*x;
      else if (idx == 33  ) num = 0.98440+0.00000*x;
      else if (idx == 34  ) num = 0.98109+0.00010*x;
      else if (idx == 35  ) num = 1.00229-0.00073*x;
      else if (idx == 36  ) num = 1.00581-0.00077*x;
      else if (idx == 37  ) num = 0.99311-0.00028*x;
      else if (idx == 38  ) num = 0.99696-0.00050*x;
      else if (idx == 39  ) num = 1.00555-0.00097*x;
      else if (idx == 40  ) num = 0.99395-0.00044*x;
      else if (idx == 41  ) num = 0.98796-0.00028*x;
      else if (idx == 42  ) num = 0.99429-0.00050*x;
      else if (idx == 43  ) num = 1.01477-0.00122*x;
      else if (idx == 44  ) num = 0.99855-0.00045*x;
      else if (idx == 45  ) num = 0.99397-0.00029*x;
      else if (idx == 46  ) num = 1.01511-0.00114*x;
      else if (idx == 47  ) num = 0.99798-0.00052*x;
      else if (idx == 48  ) num = 1.01710-0.00126*x;
      else if (idx == 49  ) num = 1.02298-0.00145*x;
      else if (idx == 50  ) num = 0.99830-0.00046*x;
      else if (idx == 51  ) num = 0.99691-0.00042*x;
      else if (idx == 52  ) num = 0.99923-0.00046*x;
      else if (idx == 53  ) num = 0.99528-0.00035*x;
      else if (idx == 54  ) num = 1.00914-0.00083*x;
      else if (idx == 55  ) num = 0.99466-0.00044*x;
      else if (idx == 56  ) num = 0.98395+0.00001*x;
      else if (idx == 57  ) num = 0.99536-0.00043*x;
      else if (idx == 58  ) num = 0.99864-0.00046*x;
      else if (idx == 59  ) num = 1.01598-0.00112*x;
      else if (idx == 60  ) num = 1.01881-0.00140*x;
      else if (idx == 61  ) num = 0.99198-0.00019*x;
      else if (idx == 62  ) num = 0.98994-0.00025*x;
      else if (idx == 63  ) num = 0.99481-0.00039*x;
      else if (idx == 64  ) num = 0.99666-0.00038*x;
      else if (idx == 65  ) num = 0.99272-0.00021*x;
      else if (idx == 66  ) num = 0.99375-0.00047*x;
      else if (idx == 67  ) num = 0.99298-0.00042*x;
      else if (idx == 68  ) num = 1.01118-0.00117*x;
      else if (idx == 69  ) num = 1.00950-0.00099*x;
      else if (idx == 70  ) num = 0.98714-0.00025*x;
      else if (idx == 71  ) num = 0.99450-0.00048*x;
      else if (idx == 72  ) num = 1.01686-0.00117*x;
      else if (idx == 73  ) num = 0.98565-0.00004*x;
      else if (idx == 74  ) num = 1.00325-0.00075*x;
      else if (idx == 75  ) num = 0.98219-0.00003*x;
      else if (idx == 76  ) num = 0.99635-0.00046*x;
      else if (idx == 77  ) num = 0.99561-0.00032*x;
      else if (idx == 78  ) num = 0.99953-0.00042*x;
      else if (idx == 79  ) num = 1.02005-0.00124*x;
      else if (idx == 80  ) num = 1.00194-0.00067*x;
      else if (idx == 81  ) num = 1.01322-0.00107*x;
      else if (idx == 82  ) num = 1.00293-0.00065*x;
      else if (idx == 83  ) num = 0.99452-0.00029*x;
      else if (idx == 84  ) num = 1.01578-0.00125*x;
      else if (idx == 85  ) num = 1.00544-0.00065*x;
      else if (idx == 86  ) num = 0.99637-0.00039*x;
      else if (idx == 87  ) num = 0.98428+0.00001*x;
      else if (idx == 88  ) num = 0.99187-0.00031*x;
      else if (idx == 89  ) num = 0.99737-0.00042*x;
      else if (idx == 90  ) num = 1.00048-0.00060*x;
      else if (idx == 91  ) num = 1.02340-0.00147*x;
      else if (idx == 92  ) num = 0.99232-0.00027*x;
      else if (idx == 93  ) num = 1.00800-0.00094*x;
      else if (idx == 94  ) num = 0.98921-0.00029*x;
      else if (idx == 95  ) num = 1.00956-0.00084*x;
      else if (idx == 96  ) num = 1.00888-0.00091*x;
      else if (idx == 97  ) num = 0.99272-0.00018*x;
      else if (idx == 98  ) num = 1.01214-0.00111*x;
      else if (idx == 99  ) num = 1.01610-0.00107*x;
      else if (idx == 100 ) num = 0.99436-0.00021*x;
      // // ADD THE FUNCTIONS FOR SYSTS BELOW! +1 SIGMA (IDX==-1) AND -1 SIGMA (IDX==-2)
      else if (idx == -1  ) num = 1.00129-0.00053*x;
      else if (idx == -2  ) num = 1.00014-0.00052*x;
   }

   if (idx==200) den = 1.;
   if (idx==300) num = den*den;

   return num / den;
}

#endif //#ifndef tnp_weight_h
