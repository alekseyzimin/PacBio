// This program finds where we shouldn't join sub-mega reads for given
// pacbio reads
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <set>
#include <algorithm>
#include <charb.hpp>
#include <misc.hpp>
#include <src_mega_reads/findGapsInCoverageOfPacbios_cmdline.hpp>

struct overlapInfoStruct {
     int impliedBegin;
     int impliedEnd;
     int actualMatchBeginOfImpliedOverlap;
     int actualMatchEndOfImpliedOverlap;
};

std::vector<struct overlapInfoStruct> overlapInfo;
std::vector<int> beginsToCover, endsToCover;
int debug;

void reportNonOverlappedGaps (int minOvl, charb &pacbio);
void createGapsToCover (int minOvlForMatchingRgns);
FILE *Fopen (const char *fn, const char *mode);
bool mySortFunction (int a, int b)
{
     if (overlapInfo[a].impliedBegin != overlapInfo[b].impliedBegin)
          return (overlapInfo[a].impliedBegin < overlapInfo[b].impliedBegin);
     return (overlapInfo[b].impliedEnd < overlapInfo[a].impliedEnd);
}

bool overlapInfoSortFunction (struct overlapInfoStruct a, struct overlapInfoStruct b)
{
     if (a.actualMatchBeginOfImpliedOverlap != b.actualMatchBeginOfImpliedOverlap)
	  return (a.actualMatchBeginOfImpliedOverlap < b.actualMatchBeginOfImpliedOverlap);
     return (a.actualMatchEndOfImpliedOverlap < b.actualMatchEndOfImpliedOverlap);
}

int main (int argc, char **argv)
{
     charb fname;
#if 0
     // The following was used for the runs up to 9/18/14
     strcpy (fname, "/genome3/raid/alekseyz/PB_ScerW303/mega-reads/coords_m300_k15_10x_70_B15_U1_mm");
     // The next was used for the 9/24 and 9/25/14 non-partial matches runs
     strcpy (fname, "mr5_1.70.13.25.0.1.blasr.out");
     // The next was used for the 9/25/14 partial matches runs
     strcpy (fname, "/home/alekseyz/mikeStuff/mrN.70.13.25.0.1.blasr.out");
     // The next was used for the 9/29/14 runs
     strcpy (fname, "/home/alekseyz/mikeStuff/mrN6_max.70.13.25.0.1.blasr.out");
     // The next was used for the 10/2/14 runs
     strcpy (fname, "test.blasr.out");
#endif
     int pacbioLen = 0;
     charb line(100), pacbio(100);
     std::set<std::string> pacbioNames;
     bool isFirstLine = true;
     FILE *infile = Fopen (fname, "r");
     fseek (infile, 0L, SEEK_END);
     long fsize = ftell (infile);
     rewind (infile);
     if (fsize == 0L) {
	  fclose (infile);
	  fprintf (stderr, "Input file %s is empty. Bye!\n", (char *) fname);
	  exit (1);
     }
     cmdline_parse args;
     args.parse (argc, argv);
     strcpy (fname, args.input_file_arg);
     int minOvlForMatchingRgns = args.max_gap_overlap_arg; // Used for the minimum actual match overlap that is not considered a gap
     int minOvl = args.min_ovl_implied_vs_gap_end_arg; // Used for the minimum overlap for the implied overlaps
     int minMatchLenForImpliedMatch = args.min_match_len_for_implied_match_arg;
     std::vector<char *> flds;
     debug = 0;
     fgets (line, 100, infile); // Skip the first line (may have to change later)
     while (fgets (line, 100, infile)) {
	  if (debug > 1) fputs (line, stdout);
	  int numFlds = getFldsFromLine (line, flds);
	  if (numFlds < 12) {
	       fprintf (stderr, "Line has %d fields, must have at least 12, line is", numFlds);
	       for (int i=0; i<numFlds; ++i) {
		    if (i > 0)
			 fprintf (stderr, " ");
		    fprintf (stderr, "%s", flds[i]);
	       }
	       fprintf (stderr, "\n");
	       exit (1);
	  }
	  if (strcmp (flds[0], (char *)pacbio) != 0) {
	       if ((pacbio.len() > 0) && (! isFirstLine)) {
		    std::sort (overlapInfo.begin(), overlapInfo.end(), overlapInfoSortFunction);
		    createGapsToCover (minOvlForMatchingRgns);
		    reportNonOverlappedGaps (minOvl, pacbio);
	       }
	       strcpy (pacbio, flds[0]);
	       std::string pacbioName = std::string (pacbio);
	       if (pacbioNames. find (pacbioName) != pacbioNames.end()) {
		    fprintf (stderr, "Pacbio read %s has records in multiple places. Bye!\n", (char *) pacbio);
		    exit (1);
	       }
	       pacbioNames.insert (pacbioName);
	       beginsToCover.clear();
	       endsToCover.clear();
	       overlapInfo.clear();
	       isFirstLine = true;
	  }
	  if (isFirstLine) {
	       pacbioLen = atoi (flds[11]);
	       isFirstLine = false; }
	  
	  int flds0 = atoi (flds[9]), flds1 = atoi (flds[10]);
	  
	  if (flds1 - flds0 >= minMatchLenForImpliedMatch) {
	       int impliedBegin = flds0 - atoi (flds[6]);
	       if (impliedBegin < 0)
		    impliedBegin = 0;
	       int impliedEnd = flds1 + (atoi (flds[8]) - atoi (flds[7]));
	       if (impliedEnd > pacbioLen)
		    impliedEnd = pacbioLen;
	       struct overlapInfoStruct oIS;
	       oIS.impliedBegin = impliedBegin;
	       oIS.impliedEnd = impliedEnd;
	       oIS.actualMatchBeginOfImpliedOverlap = flds0;
	       oIS.actualMatchEndOfImpliedOverlap = flds1;
	       overlapInfo.push_back (oIS);
	  }
     }

     if ((pacbio.len() > 0) && (! isFirstLine)) {
	  std::sort (overlapInfo.begin(), overlapInfo.end(), overlapInfoSortFunction);
	  createGapsToCover (minOvlForMatchingRgns);
	  reportNonOverlappedGaps (minOvl, pacbio);
     }

     return (0);
}

void reportNonOverlappedGaps (int minOvl, charb &pacbio)
{
     int i, j, k;
     int spclOvlVal = 1;
     std::vector<int> indices, reducedImpliedBegins, reducedImpliedEnds;
     std::vector<int> emptyIntVector;
     emptyIntVector.clear();
     if (overlapInfo.size() == 0)
	  return;
     if (beginsToCover.size() == 0)
	  return;
     for (i=0; i<(int) overlapInfo.size(); ++i)
	  indices.push_back (i);
     std::sort (indices.begin(), indices.end(), mySortFunction);
     std::vector<std::vector<int> >vectorOfMatchRgnsThatKillGapRgn;
     vectorOfMatchRgnsThatKillGapRgn.clear();
     for (i=0; i<(int) beginsToCover.size(); ++i)
	  vectorOfMatchRgnsThatKillGapRgn.push_back (emptyIntVector);
     for (i=0; i<(int) indices.size(); ++i) {
	  if (debug > 1) std::cout << "i = " << i << "\n";
	  if (overlapInfo[indices[i]].impliedEnd - overlapInfo[indices[i]].impliedBegin <= 2 * minOvl)
	       continue;
	  for (j=0; j<(int) beginsToCover.size(); ++j) {
	       if (debug > 1) std::cout << "j = " << j << " rgn = (" << beginsToCover[j] << "," << endsToCover[j] << ")\n";
	       if (overlapInfo[indices[i]].impliedBegin > beginsToCover[j] - minOvl)
		    continue;
	       if (overlapInfo[indices[i]].impliedEnd < endsToCover[j] + minOvl)
		    break;
	       // We are covering the current gap with the implied ovl
	       bool priorFound=false, followFound=false;
	       for (k=j; k>0; --k) {
		    if ((overlapInfo[indices[i]].actualMatchBeginOfImpliedOverlap <= beginsToCover[k] - spclOvlVal) && (overlapInfo[indices[i]].actualMatchEndOfImpliedOverlap >= endsToCover[k-1] + spclOvlVal)) {
			 vectorOfMatchRgnsThatKillGapRgn[j].push_back (k);
			 if (debug > 1) std::cout << "Rgn to kill at 1: (" << beginsToCover[j] << "," << endsToCover[j] << "); ";
			 if (debug > 1) std::cout << "Killing rgn: (" << overlapInfo[indices[i]].actualMatchBeginOfImpliedOverlap << "," << overlapInfo[indices[i]].actualMatchEndOfImpliedOverlap << ") in (" << overlapInfo[indices[i]].impliedBegin << "," << overlapInfo[indices[i]].impliedEnd << ")\n";
			 priorFound = true;
			 break; } }
	       k=0;
	       if ((! priorFound) && (overlapInfo[indices[i]].actualMatchBeginOfImpliedOverlap <= beginsToCover[k] - spclOvlVal)) {
		    vectorOfMatchRgnsThatKillGapRgn[j].push_back (0);
		    if (debug > 1) std::cout << "Rgn to kill at 1.2: (" << beginsToCover[j] << "," << endsToCover[j] << "); ";
		    if (debug > 1) std::cout << "Killing rgn: (" << overlapInfo[indices[i]].actualMatchBeginOfImpliedOverlap << "," << overlapInfo[indices[i]].actualMatchEndOfImpliedOverlap << ") in (" << overlapInfo[indices[i]].impliedBegin << "," << overlapInfo[indices[i]].impliedEnd << ")\n";
		    priorFound = true;
	       }
	       for (k=j+1; k<(int) beginsToCover.size(); ++k) {
		    if (debug > 1) std::cout << "k = " << k << "\n";
		    if ((overlapInfo[indices[i]].actualMatchBeginOfImpliedOverlap <= beginsToCover[k] - spclOvlVal) && (overlapInfo[indices[i]].actualMatchEndOfImpliedOverlap >= endsToCover[k-1] + spclOvlVal)) {
			 vectorOfMatchRgnsThatKillGapRgn[j].push_back (k);
			 if (debug > 1) std::cout << "Rgn to kill at 2: (" << beginsToCover[j] << "," << endsToCover[j] << "); ";
			 if (debug > 1) std::cout << "Killing rgn: (" << overlapInfo[indices[i]].actualMatchBeginOfImpliedOverlap << "," << overlapInfo[indices[i]].actualMatchEndOfImpliedOverlap << ") in (" << overlapInfo[indices[i]].impliedBegin << "," << overlapInfo[indices[i]].impliedEnd << ")\n";
			 followFound = true;
			 break; } }
	       k=(int) beginsToCover.size()-1;
	       if ((! followFound) && (overlapInfo[indices[i]].actualMatchEndOfImpliedOverlap >= endsToCover[k] + spclOvlVal)) {
		    vectorOfMatchRgnsThatKillGapRgn[j].push_back (k+1);
		    if (debug > 1) std::cout << "Rgn to kill at 2.2: (" << beginsToCover[j] << "," << endsToCover[j] << "); ";
		    if (debug > 1) std::cout << "Killing rgn: (" << overlapInfo[indices[i]].actualMatchBeginOfImpliedOverlap << "," << overlapInfo[indices[i]].actualMatchEndOfImpliedOverlap << ") in (" << overlapInfo[indices[i]].impliedBegin << "," << overlapInfo[indices[i]].impliedEnd << ")\n";
		    followFound = true;
	       }
	  }
     }
     int intervalBegin = -1, intervalEnd = -1;
     for (i=0; i<(int) beginsToCover.size(); ++i) {
	  if ((int) vectorOfMatchRgnsThatKillGapRgn[i].size() > 1) {
	       std::sort (vectorOfMatchRgnsThatKillGapRgn[i].begin(), vectorOfMatchRgnsThatKillGapRgn[i].end());
	       if (debug) std::cout << "Vector of match rgns that kill region " << i << ": ";
	       for (j=0; j<(int) vectorOfMatchRgnsThatKillGapRgn[i].size(); ++j)
		    if (debug) std::cout << " " << vectorOfMatchRgnsThatKillGapRgn[i][j];
	       if ((vectorOfMatchRgnsThatKillGapRgn[i][0] <= i) && (vectorOfMatchRgnsThatKillGapRgn[i][vectorOfMatchRgnsThatKillGapRgn[i].size()-1] > i)) {
		    if (debug) std::cout << "  IT'S DEAD\n";
		    if (intervalBegin < 0) {
			 intervalBegin = beginsToCover[i];
			 intervalEnd = endsToCover[i]; }
		    if (beginsToCover[i] > intervalEnd) {
			 std::cout << pacbio << " " << intervalBegin << " " << intervalEnd << "\n";
			 intervalBegin = beginsToCover[i]; }
		    if (endsToCover[i] > intervalEnd)
			 intervalEnd = endsToCover[i];
	       }
	       if (debug) std::cout << "\n";
	  }
     }
     if (intervalBegin >= 0)
	  std::cout << pacbio << " " << intervalBegin << " " << intervalEnd << "\n";
}

void createGapsToCover (int minOvlForMatchingRgns)
{
     int begin = 0, end = 0;
     beginsToCover.clear();
     endsToCover.clear();
     for (int i=0; i<(int) overlapInfo.size(); ++i) {
	  begin = overlapInfo[i].actualMatchBeginOfImpliedOverlap;
	  if (begin > end - minOvlForMatchingRgns) {
	       if (end > 0) {
		    int firstToCover, lastToCover;
		    if (end < begin) {
			 firstToCover = end;
			 lastToCover = begin; }
		    else {
			 firstToCover = begin;
			 lastToCover = end; }
		    int vSize = beginsToCover.size();
		    if ((vSize == 0) || (firstToCover != beginsToCover[vSize-1]) || (lastToCover != endsToCover[vSize-1])) {
			 beginsToCover.push_back (firstToCover);
			 endsToCover.push_back (lastToCover); }
	       }
	  }
	  if (end < overlapInfo[i].actualMatchEndOfImpliedOverlap)
	       end = overlapInfo[i].actualMatchEndOfImpliedOverlap;
     }
}

FILE *Fopen (const char *fn, const char *mode)
{
     FILE *result;
     result = fopen (fn, mode);
     if (result == NULL)
     {
          fprintf (stderr, "Couldn't open file '%s' for ", fn);
          switch (mode[0])
          {
          case 'r':
               fprintf (stderr, "reading");
               break;
          case 'w':
               fprintf (stderr, "writing");
               break;
          case 'a':
               fprintf (stderr, "appending");
               break;
          default:
               fprintf (stderr, "unknown operation code '%c'", mode[0]);
               break;
          }
          fprintf (stderr, ". Bye!\n");
          exit (-1);
     }

     return (result);
}

