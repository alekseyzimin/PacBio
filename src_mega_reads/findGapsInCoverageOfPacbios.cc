// This program find where we might be able to join sub-mega reads for given
// pacbio reads
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <algorithm>
#include <charb.hpp>
#include <misc.hpp>

std::vector<int> impliedBegins, impliedEnds, beginsToCover, endsToCover;

void reportNonOverlappedGaps (int minOvl, charb &pacbio);
FILE *Fopen (const char *fn, const char *mode);
bool mySortFunction (int a, int b)
{
     if (impliedBegins[a] != impliedBegins[b])
	  return (impliedBegins[a] < impliedBegins[b]);
     return (impliedEnds[b] < impliedEnds[a]);
}

int main (int argc, char **argv)
{
     charb fname;
     strcpy (fname, "/genome3/raid/alekseyz/PB_ScerW303/mega-reads/coords_m300_k15_10x_70_B15_U1_mm");
     int minOvlForMatchingRgns = 100;
     int minOvl = 0;
     int begin = 0, end = 0, pacbioLen = 0;
     int minMatchLenForImpliedMatch = 25;
     charb line(100), pacbio(100);
     bool isFirstLine = true;
     FILE *infile = Fopen (fname, "r");
     std::vector<char *> flds;
     while (fgets (line, 100, infile)) {
	  int numFlds = getFldsFromLine (line, flds);
	  if ((numFlds > 0) && (flds[0][0] == '>')) {
	       if ((pacbio.len() > 0) && (! isFirstLine))
		    reportNonOverlappedGaps (minOvl, pacbio);
	       strcpy (pacbio, flds[1]);
	       begin = end = 0;
	       beginsToCover.clear();
	       endsToCover.clear();
	       impliedBegins.clear();
	       impliedEnds.clear();
	       isFirstLine = true;
	       continue; }
	  if (isFirstLine) {
	       pacbioLen = atoi (flds[9]);
	       isFirstLine = false; }
	  
	  int flds0 = atoi (flds[0]), flds1 = atoi (flds[1]);
	  
	  begin = flds0 - 1;
	  if (flds1 - flds0 + 1 >= minMatchLenForImpliedMatch) {
	       int impliedBegin = flds0 - atoi (flds[2]);
	       if (impliedBegin < 0)
		    impliedBegin = 0;
	       int impliedEnd = flds1 + (atoi (flds[10]) - atoi (flds[3]));
	       if (impliedEnd > pacbioLen)
		    impliedEnd = pacbioLen;
	       impliedBegins.push_back (impliedBegin);
	       impliedEnds.push_back (impliedEnd);
	  }
	  if (begin > end - minOvlForMatchingRgns) {
	       if (end > 0) {
		    if (end < begin) {
			 beginsToCover.push_back (end);
			 endsToCover.push_back (begin); }
		    else {
			 beginsToCover.push_back (begin);
			 endsToCover.push_back (end); }
	       }
	  }
	  if (end < flds1)
	       end = flds1;
     }

     if ((pacbio.len() > 0) && (! isFirstLine))
	  reportNonOverlappedGaps (minOvl, pacbio);

     return (0);
}

void reportNonOverlappedGaps (int minOvl, charb &pacbio)
{
     int i, j, intervalBegin, intervalEnd;
     std::vector<int> indices, reducedImpliedBegins, reducedImpliedEnds;
     if (impliedBegins.size() == 0)
	  return;
     if (beginsToCover.size() == 0)
	  return;
     for (i=0; i<impliedBegins.size(); ++i)
	  indices.push_back (i);
     std::sort (indices.begin(), indices.end(), mySortFunction);
     // spcl => mySortFunction; check above syntax
     reducedImpliedBegins.push_back (impliedBegins[indices[0]]);
     reducedImpliedEnds.push_back (impliedEnds[indices[0]]);
     for (i=1; i<impliedBegins.size(); ++i) {
	  if (impliedEnds[indices[i]] <= reducedImpliedEnds[reducedImpliedEnds.size()-1])
	       continue;
	  reducedImpliedBegins.push_back (impliedBegins[indices[i]]);
	  reducedImpliedEnds.push_back (impliedEnds[indices[i]]);
     }
     i = j = 0;
     intervalBegin = intervalEnd = -1;
     while (1) {
	  if (reducedImpliedEnds[j] < endsToCover[i] + minOvl) {
	       ++j;
	       if (j >= reducedImpliedEnds.size())
		    break;
	       continue; }
	  if (reducedImpliedBegins[j] <= beginsToCover[i] - minOvl) {
	       if (intervalBegin < 0) {
		    intervalBegin = beginsToCover[i];
		    intervalEnd = endsToCover[i]; }
	       if (beginsToCover[i] > intervalEnd) {
		    std::cout << "BREAK " << pacbio << " " << intervalBegin << " " << intervalEnd << "\n";
		    intervalBegin = beginsToCover[i]; }
	       if (endsToCover[i] > intervalEnd)
		    intervalEnd = endsToCover[i];
	  }
	  ++i;
	  if (i >= beginsToCover.size())
	       break;
     }
     if (intervalBegin >= 0)
	  std::cout << "BREAK " << pacbio << " " << intervalBegin << " " << intervalEnd << "\n";
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

