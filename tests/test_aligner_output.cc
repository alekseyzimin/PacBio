#include <gtest/gtest.h>
#include <src_jf_aligner/superread_parser.hpp>
#include <src_jf_aligner/pb_aligner.hpp>

namespace {
struct remove_file {
  const char* path;
  bool do_unlink;
  remove_file(const char* p, bool unlink = true) : path(p), do_unlink(unlink) { }
  ~remove_file() { if(do_unlink) unlink(path); }
};

static const char* const pac_bio =
  ">pb\n"
  "CGCGAAAGGATCGCCTCGATCGCTGGAGCCAAGCGCGTGGCCGGTCCGAGGGCGGTCGCGGTAGGCGTCG"
  "GCCAGATCGTCGGGCCAGGTTCCGGCGGCGCCAGCAGACTCGGTCCAGAGGAAGGGCAGCCCCTGCACCG"
  "CGATTCCTTCCCGCGCCTTGGCCGGCCGGCCGTCGCGGCAGGCGGTACGGCCACGAATCGCTGTGGATTG"
  "CGAACCGCATCATCTCTCTATCAACAACAACAACGGAAGGAGGAAGGAAAAGGAGAGAGATGAGCGGTTC"
  "GCGAATCCACAGCGATCGTGGCGACCGCCGCCGCGACGGCCCAAAGGACCGGGCCAAGGCGGCGAAGGAT"
  "CGCGGGCAGGCGGCTTCCTGAACCGAGCTGCTGGCCGCCGACACGCCCGGAACGATCTGGCCGACGCCTA"
  "ACACGTCGACGCCCCTCCCGGAACGCGGGCCAAGCGCTGGCGCCAGCGGATCGAGGCGAGCCTTTCGCGG"
  "ATGCCCAGCCACGGCGGAGATTCTTGCGCGACTCGGCTGAGAGAGACCGCTCCGCGCCGGGGCCCCTATG"
  "CCTTGCCGCGGGGTTCGTCGGGAGGAAGGGGAACTTGGCAGCGGGACGCCGCGCTGCCGGCAAAGGCGAG"
  "CGCACCCCCCGCGGAGAGGCGGAGGTCCTTGCCGATCCCGCCGGGGGCAGGGAGGCGGCTGCGGCCGGGA"
  "AGAGGGGTCATCCAACCCCGGGAGAGGCCGGATGTGCCCTGCCGGATCCGGGCAGGGCGCCGCCCGGGTC"
  "CCAGTAGGCCCCGGCGAAGAAGAGCGCAACCATATCGGCCAGAGGCATCGGGACCGGCAGGATTCGGTGT"
  "AGGAGGCGGCCGGCGTCCGCCCGACCAGGCATGGCCGATCCCCGATGCGCCCAGAGTTCGGGCCAGTGCG"
  "GTTTCCGTCGGCAGCCGAAGGGTCCTCGCGTTCCGCCATTTGTCCGCCCCGCCGGGGCTGCGCCCCTGCT"
  "CGTTGGCCGGCGCCGCACGTGCGCGCCCACCGGACGCCGCGCGACCAGCCCGCTTGCGCAATGTCGGGAT"
  "GCACCGTGGCAATCCGCGTGCCGTGGAAAGACGAGGTGCGCGGCGCCCGCCCTGCCCGCGGGCCGTTGCC"
  "GCCTGCCGTTCGCCCCGCATGGCCGCGAAAGGCCGAGATTGACATCTACAAAAACAAAAAATGCCTCGCC"
  "GACCGGTGGCGTCGCCCTTCACCATACGCACGCAAACCAACGCACACCACCCAAGCCTGTCTCTACCCTC"
  "CGCCACCCACCCACCGCAGATCGTACC\n";
static const char* const super_reads =
  ">234R_239F\n"
  "GAGCCGCCGGCGTAGATGATCGTCTCGACATTGGCGACCTCGCCCAGATCGAGCCCCTCGGTCAGCACGA"
  "TGATCGTGTCGTTCCCGCCGTCGGCGAATTCGCGCACGATATCGTCGGCGGAATCGACGAAATAGGTGTC"
  "GTCGCCGCCGAAGCCCGCCATCCGGTCGCCGCCGCCGCGCCCGTCGAGCACGTTGTCGCCGTCATTGCCG"
  "ATGAGCGTATCGGCGAAGCGGCTGCCGATGCCCGCCTCGATGTAATAGTCGCTGGCGTCCTCGTTGTGGC"
  "CCGCGAAGATCGAGACGTTGTTGCTGTGGCCGTTCACGCTCGAGAACCGGCCTGCGCGCAGATCGAGGAT"
  "CGTGCCCGCCGTCGAGCCCGAGAAGTCGATCGTGTCCATGCCGCCGGTGTCGTGGATCACGAAGCCGATG"
  "TCGTGCTGCAGCCCCTCGGACATCAGCTCGTAGCCGAAGACGCTGCTGTTGAAGCCGTAGATGTCGTCGC"
  "CCGTGTTGAGGTCGATGGTCTGGTAGGTCCGCACCCCGTCCTCGTCCACCGTCGAGAAGAAGCGGCGGAT"
  "CACGGCCTCGATATCGGCCATGAGCGGCGTGGTGGCGGACCAGCGGCTCGTCTCGCCCACGTCGATGCCG"
  "TCGAAATAGGACATGACGCTGTATTGCTGCCGGTCGTAGGTCCAGGTCGCGTTGTTCAGATAGTTGATCT"
  "GCACGCCGCCGGGACCGCTGTAATTGTAGAGGCCGGGGTGGTTCAGGCCGAATTCATGGCCGAATTCATG"
  "GACGAACGAGTCGAACACATAGCCGCCGAGAGCGGTCTTGTCGGGCTCGGTGTCGTGAAAGCGCTGGCCG"
  "ATGCTGACATAGCGGTTGCTCGAGAAGGCGCTGCCGTCGCTCGGTTCGCCGATCTCGGGGCTCACCACCT"
  "CCATCCAGTCCGTCGCGCGGTCGAAGGGCGCGTCGTCCACGATCTCGAACTGGAGCGGCGTGACCGTGGC"
  "CCAGGTCTGCAGCGCGCCGATGGCGGCCGCCTTGTAATCCGGGTGGTCGTTCAGGCCGGTCACGTTGATC"
  "CGGAGCGGCTCGGCCGCGATGGCCGGGGTGAGGCTGCGGTTCAGCGCCCCCAGATAATCGAGGAAGGCAT"
  "CGTAGTCCTCGCTCTTGCCATAAGGATTATCGGGCGTCTGATCGTTCACCCACCATGTGTCCGCGACTTG"
  "TCTGGCAGACGGCATCTTCGGCTCTCCCTGACTGAAGGCTGGCTAGGCCCATGCTTCGGTGCTGCCGAAC"
  "AGCGTTGCCCCGTTCCTGTTCCCGTTCAAGGAAAGTCAAACCGGCGTGATCTTGCCTATGGGTTGAAAAA"
  "CCGCCCCGTGACGTGCCCGGAATTTCATGAAGAATGAAAGCTCAAGCCCTGAGTGCTTCCCCGCCGGCGA"
  "GCTTCATCAGCAGCAGGATGACGTTGCGTAGCCACGGTTTCCCGCTTTCGTCTCACGCTGGCGCGACGGC"
  "GGGACGGCGGGTCTTCCGGCCGCGGACGCGCCGGTTCGGGCCGCCGACCAGACGGACCCTCCGCCGCTTC"
  "GACTGTATGGGCCGCAGCGTCGGAACAGAGGCGGCGGGGCAGGGTTGCCCGGGCACGACGCCGCATTCGA"
  "AAGGAGGACCGATGAGCGACCCCCTGATCTTCGACAGCTGGACCGGCCTCGGGCGGATCCTGCTGGCCGG"
  "CACCCTCGCCTACGGCGCGTTGGTGGCGATGCTGCGCGTCTCGGGCAAGCGCACGCTGAGCAAGATGAAC"
  "GCCTTCGACCTGATCGTGACCGTGGCCCTCGGCTCGACGCTGGCCACCGTGATCCTGAACCGGAGCGTGC"
  "CGCTGGCCGAGGGCGTGCTGGCGCTGGCCCTGCTGATCTGCCTGCAATATGCGATCACCTGGACATCGGT"
  "GCGCTGGCCTCTGGCCGAGAGCATCGTCAAGAGCGAGCCCACGCTGCTCCTGCACGAAGGCTGCTTCCTC"
  "GACGATGCCATGCGCCGTCAGCGCGTCACGCAGGCCGAGATCCTCGCGGCCCTGCGCGAGAGCGGCCTCG"
  "ACAGCCCCGGCGCCGCCCGCTCCGTGGTGCTCGAGACGGACGGCAGCCTCTCGGTCGTGTCGCGCACGGA"
  "CTGACGCGGCGTGGAGGCGGGCTCTGGGCCGTGTCCGCCGTGTGCATCAAGAAAATTCCCGCTGCACACT"
  "TGCTCGATGTGCAGCGGAAAGATATTTGATGCACACACTGGCCCTGAGGCGCAGACGGTCGGAGGAGACC"
  "CGCCATGACGATCCGCTCGCATTCCCGCGCCGCGCAGGTCGCCTATCAGGACCTGCTCCGGCTTCATCTC"
  "GACGAGAGCGCCTCGGCGCTGATCGGCAGCATCGAGCAGCGCACCCGCAACGGGCGGATCTATCTCTACG"
  "ACAAGTTCCGCATCGGCACCGAGATGAAGAGCCGCTACCTCGGCGAGGGCACGCCCGAGCTTGTCGCGCG"
  "GCTCGAGCGGGCGGCGGCGCTGAAGGCCGGGGCCGAGGAGCGCCGGGCGACCATGGCGCGCTTGGCCCGG"
  "GTGCTGCGCGCGGAGGGCTTCACCGGCACCGACCGCGAGACGGGCTCGCTGCTTCTGGCCTTCGCGCGGG"
  "CGGGCGTGTTCCGGCTGGGCGGGACGCTGGTGGGCACGGCCGCCTATGCGCTCTATCAGGGCGAGCTCGG"
  "CGTGCGGTTCGATGCCGAGGAACTGGCCCAGACCGGGGACATCGACTTCGCGAGTTTCGAGCGGCTGTCG"
  "GTCGCGCTCGGCGACCGCGTCGAGGAGGAGCCGGGCGACATCCTGCAGGCGCTGAAGTTCGATCCGGTGC"
  "CGGGGTTGGCCGACCGGCAGGTCTGGAAATGGCGCCAGAGCCGCGGGCAGGCGATGGTCGAGTTCCTGAC"
  "GCCCGCCTTCGGCGACGAACGGGTGAAGCCGCTGCCCGCGCTCGGCATCAGTGCGCAGGCGCTGAACTAT"
  "CTCAACTTCCTGATCGCGGAGCCCATTCCCGCCGTGGCGCTCTACCGCTCCGGCGTGCTGGTGCAGATCC"
  "CGCGGCCCGAGCGGTTCGCGATCCACAAGCTGATCGTGGCCGACCGCCGCCGCGACGGCCCGGACCGGGC"
  "CAAGGCGCGGAAGGATCGCGGGCAGGCGGCCTTCCTGACCGAGCTGCTGGCCGCCGACCGGCCCGACGAT"
  "CTGGCCGACGCCTACCGCGACGCCCTCGGACGCGGCCCGCGCTGGCGCCAGCGGATCGAGGCGAGCCTTT"
  "CGCGGATGCCCGCCACGGCGGAGATCCTGCGCGACCTCGGCTGAGAGAGACCGCTCCCGCGCCGGGCCTA"
  "TGCCTGCCGCGGGGTCTCGTGCGGGAGGAGGGGGCCTTGGCGCGGGACGCCGCCGCTGCCGGCAGGCGAG"
  "GCGGCACCCCGCGGAGAGGCGGAGGTCCTTGCCGATCCCGCCGGGGCAGGGAGGCGGCTGGCGGCGGGCA"
  "AGAGGGCATCCCCCCGGGAGAGGCGGATGCCCTGCCGATCCGGGCAGGGCGCCGCCCCGGGTCCGCTAGG"
  "CCCCGGCGAGGAAGAAGCGCACCATCTCGGCCGAGGCATCGGGACCGGCCGGATCGGTGTAGGAGCCGGC"
  "GGCCCGTCCGCCCGACCAGGCATGGCCCGATCCCTCGATGCGCCAGAGTTCGGCCAGTGCGGTTCCGTCG"
  "GCAGCCGCATGGGTCTCGCGTCGCCATG\n"
  ">98345F_1234567R_5472F\n"
  "GCGCAGGCCGCGCCGGAGGCCTCCGCCGCCCCCCCGCGTCGCCGGATGGGGGCCGCGGTCGATGCCCTCC"
  "GCCGGCTGCGCCCCACGGGCCTCGCGTCCGCCGACGCCGGACGGGGGGCCGAGCCTGCCCTTCCCACGGG"
  "CGCCCGGTTCGAGCGGCTCCATCACGCAGGCACCGCCGGTGCCCGCGACTACCGGCTCTATATACCCGCG"
  "GCGCAGCCGGAGGGTCGGCCGGGGCTCGTCCTGATGCTGCACGGCTGCACCCAGACGCCCGAGGATTTCG"
  "CGGCAGGCACCGGGATGAACGCGGCGGCCGAGCGGCACGGGCTGATCGTCCTCTATCCGCATCAGGACCG"
  "CAGCCACAATGCGCAGGGCTGCTGGAACTGGTTCCGGCCCGGCGATCAGGAGGCCCACCGCGGCGAGCCC"
  "GCGCTTCTGGCCGACCTCGTGCGGGCGGTGGCCGCCCGGCATGACGTGGATGCGGGGCGGATCTTCGTCG"
  "CGGGCCTCTCGGCGGGCGGGGCGATGGCCGCGACGCTCGGGGCCACGCACCCCGAGCTCTTCTCGGCCGT"
  "GGGGGTCCATTCCGGCCTGCCGCACCGGGCGGCGCACGATGTCATCTCGGCCTTCGCGGCCATGCGGGGC"
  "GACGGGCAGGCGGCACGGCCCGCGGCAGGCGGGGCGCCGCGCACCATCGTCTTCCACGGCACGGCGGATG"
  "CCACGGTGCATCCCGACAATGCGCAGCGGCTGATCGACGCGGCGGTCGGCGGGCGGCAGTCGGCGGCGCG"
  "CACCGAGCAGGGGCGCAGCCCCGGCGGGCGGACATGGCGACGCGAGACCCATGCGGCTGCCGACGGAACC"
  "GCACTGGCCGAACTCTGGCGCATCGAGGGATCGGGCCATGCCTGGTCGGGCGGACGGGCCGCCGGCTCCT"
  "ACACCGATCCGGCCGGTCCCGATGCCTCGGCCGAGATGGTGCGCTTCTTCCTCGCCGGGGCCTAGCGGAC"
  "CCGGGGCGGCGCCCTGCCCGGATCGGCAGGGCATCCGCCTCTCCCGGGGGGATGCCCTCTTGCCCGCCGC"
  "CAGCCGCCTCCCTGCCCCGGCGGGATCGGCAAGGACCTCCGCCTCTCCGCGGGGTGCCGCCTCGCCTGCC"
  "GGCAGCGGCGGCGTCCCGCGCCAAGGCCCCCTCCTCCCGCACGAGACCCCGCGGCAGGCATAGGCCCGGC"
  "GCGGGAGCGGTCTCTCTCAGCCGAGGTCGCGCAGGATCTCCGCCGTGGCGGGCATCCGCGAAAGGCTCGC"
  "CTCGATCCGCTGGCGCCAGCGCGGGCCGCGTCCGAGGGCGTCGCGGTAGGCGTCGGCCAGATCGTCGGGC"
  "CGGTCGGCGGCCAGCAGCTCGGTCAGGAAGGCCGCCTGCCCGCGATCCTTCCGCGCCTTGGCCCGGTCCG"
  "GGCCGTCGCGGCGGCGGTCGGCCACGATCAGCTTGTGGATCGCGAACCGCTCGGGCCGCGGGATCTGCAC"
  "CAGCACGCCGGAGCGGTAGAGCGCCACGGCGGGAATGGGCTCCGCGATCAGGAAGTTGAGATAGTTCAGC"
  "GCCTGCGCACTGATGCCGAGCGCGGGCAGCGGCTTCACCCGTTCGTCGCCGAAGGCGGGCGTCAGGAACT"
  "CGACCATCGCCTGCCCGCGGCTCTGGCGCCATTTCCAGACCTGCCGGTCGGCCAACCCCGGCACCGGATC"
  "GAACTTCAGCGCCTGCAGGATGTCGCCCGGCTCCTCCTCGACGCGGTCGCCGAGCGCGACCGACAGCCGC"
  "TCGAAACTCGCGAAGTCGATGTCCCCGGTCTGGGCCAGTTCCTCGGCATCGAACCGCACGCCGAGCTCGC"
  "CCTGATAGAGCGCATAGGCGGCCGTGCCCACCAGCGTCCCGCCCAGCCGGAACACGCCCGCCCGCGCGAA"
  "GGCCAGAAGCAGCGAGCCCGTCTCGCGGTCGGTGCCGGTGAAGCCCTCCGCGCGCAGCACCCGGGCCAAG"
  "CGCGCCATGGTCGCCCGGCGCTCCTCGGCCCCGGCCTTCAGCGCCGCCGCCCGCTCGAGCCGCGCGACAA"
  "GCTCGGGCGTGCCCTCGCCGAGGTAGCGGCTCTTCATCTCGGTGCCGATGCGGAACTTGTCGTAGAGATA"
  "GATCCGCCCGTTGCGGGTGCGCTGCTCGATGCTGCCGATCAGCGCCGAGGCGCTCTCGTCGAGATGAAGC"
  "CGGAGCAGGTCCTGATAGGCGACCTGCGCGGCGCGGGAATGCGAGCGGATCGTCATGGCGGGTCTCCTCC"
  "GACCGTCTGCGCCTCAGGGCCAGTGTGTGCATCAAATATCTTTCCGCTGCACATCGAGCAAGTGTGCAGC"
  "GGGAATTTTCTTGATGCACACGGCGGACACGGCCCAGAGCCCGCCTCCACGCCGCGTCAGTCCGTGCGCG"
  "ACACGACCGAGAGGCTGCCGTCCGTCTCGAGCACCACGGAGCGGGCGGCGCCGGGGCTGTCGAGGCCGCT"
  "CTCGCGCAGGGCCGCGAGGATCTCGGCCTGCGTGACGCGCTGACGGCGCATGGCATCGTCGAGGAAGCAG"
  "CCTTCGTGCAGGAGCAGCGTGGGCTCGCTCTTGACGATGCTCTCGGCCAGAGGCCAGCGCACCGATGTCC"
  "AGGTGATCGCATATTGCAGGCAGATCAGCAGGGCCAGCGCCAGCACGCCCTCGGCCAGCGGCACGCTCCG"
  "GTTCAGGATCACGGTGGCCAGCGTCGAGCCGAGGGCCACGGTCACGATCAGGTCGAAGGCGTTCATCTTG"
  "CTCAGCGTGCGCTTGCCCGAGACGCGCAGCATCGCCACCAACGCGCCGTAGGCGAGGGTGCCGGCCAGCA"
  "GGATCCGCCCGAGGCCGGTCCAGCTGTCGAAGATCAGGGGGTCGCTCATCGGTCCTCCTTTCGAATGCGG"
  "CGTCGTGCCCGGGCAACCCTGCCCCGCCGCCTCTGTTCCGACGCTGCGGCCCATACAGTCGAAGCGGCGG"
  "AGGGTCCGTCTGGTCGGCGGCCCGAACCGGCGCGTCCGCGGCCGGAAGACCCGCCGTCCC\n";

static const char* normal_coords[3] = {
  "Rstart Rend Qstart Qend Nmers Rcons Qcons Rcover Qcover Rlen Qlen Rname Qname",
  "303 877 3051 3597 90 79 79 263 262 1287 3668 pb 234R_239F",
  "303 1150 1420 613 107 93 93 328 327 1287 3000 pb 98345F_1234567R_5472F"
};

static const char* forward_coords[3] = {
  "Rstart Rend Qstart Qend Nmers Rcons Qcons Rcover Qcover Rlen Qlen Rname Qname",
  "303 877 3051 3597 90 79 79 263 262 1287 3668 pb 234R_239F",
  "303 1150 1581 2388 107 93 93 328 327 1287 3000 pb 5472R_1234567F_98345R"
};

static const char* normal_details[2] = {
  "pb 234R_239F 56:-3155 57:-3154 58:-3153 59:-3152 60:-3151 61:-3150 62:-3149 63:-3148 64:-3147 65:-3146"
  " 66:-3145 67:-3144 68:-3143 69:-3142 70:-3141 71:-3140 [303:3051] [304:3052] [305:3053] [306:3054]"
  " [307:3055] [345:3090] [346:3091] [347:3092] [348:3093] [349:3094] [373:3118] [374:3119] [375:3120]"
  " [376:3121] [377:3122] [403:3146] [404:3147] [405:3148] [455:3189] [456:3190] [457:3191] [458:3192]"
  " [459:3193] [460:3194] [461:3195] [462:3196] [463:3197] [464:3198] [465:3199] [466:3200] [467:3201]"
  " [468:3202] [469:3203] [470:3204] [471:3205] [472:3206] [473:3207] [474:3208] [475:3209] [476:3210]"
  " [477:3211] [478:3212] [479:3213] [480:3214] [481:3215] [523:3256] [524:3257] [525:3258] [526:3259]"
  " [527:3260] [528:3261] [638:3367] [639:3368] [640:3369] [641:3370] [642:3371] [643:3372] [644:3373]"
  " [645:3374] [646:3375] [647:3376] [648:3377] [649:3378] [650:3379] [651:3380] [652:3381] [653:3382]"
  " [654:3383] [655:3384] [656:3385] [657:3386] [658:3387] [659:3388] [660:3389] [661:3390] [675:3403]"
  " [676:3404] [746:3467] [747:3468] [748:3469] [749:3470] [750:3471] [813:3533] [855:3575] [856:3576]"
  " [857:3577] [858:3578] [859:3579] [860:3580] [861:3581]",
  "pb 98345F_1234567R_5472F 56:1300 57:1301 58:1302 59:1303 60:1304 61:1305 62:1306 63:1307 64:1308 65:1309 66:1310 67:1311"
  " 68:1312 69:1313 70:1314 71:1315 [303:-1404] [304:-1403] [305:-1402] [306:-1401] [307:-1400] [345:-1365]"
  " [346:-1364] [347:-1363] [348:-1362] [349:-1361] [373:-1337] [374:-1336] [375:-1335] [376:-1334]"
  " [377:-1333] [403:-1309] [404:-1308] [405:-1307] [455:-1266] [456:-1265] [457:-1264] [458:-1263]"
  " [459:-1262] [460:-1261] [461:-1260] [462:-1259] [463:-1258] [464:-1257] [465:-1256] [466:-1255]"
  " [467:-1254] [468:-1253] [469:-1252] [470:-1251] [471:-1250] [472:-1249] [473:-1248] [474:-1247]"
  " [475:-1246] [476:-1245] [477:-1244] [478:-1243] [479:-1242] [480:-1241] [481:-1240] [523:-1199]"
  " [524:-1198] [525:-1197] [526:-1196] [527:-1195] [528:-1194] [638:-1088] [639:-1087] [640:-1086]"
  " [641:-1085] [642:-1084] [643:-1083] [644:-1082] [645:-1081] [646:-1080] [647:-1079] [648:-1078]"
  " [649:-1077] [650:-1076] [651:-1075] [652:-1074] [653:-1073] [654:-1072] [655:-1071] [656:-1070]"
  " [657:-1069] [658:-1068] [659:-1067] [660:-1066] [661:-1065] [675:-1052] [676:-1051] [746:-988]"
  " [747:-987] [748:-986] [749:-985] [750:-984] [813:-922] [855:-880] [856:-879] [857:-878] [858:-877]"
  " [859:-876] [860:-875] [861:-874] [959:-782] [960:-781] [961:-780] [962:-779] [963:-778] [964:-777]"
  " [965:-776] [966:-775] [967:-774] [1043:-702] [1044:-701] [1045:-700] [1046:-699] [1131:-616] [1132:-615]"
  " [1133:-614] [1134:-613]"
};

static const char* comp_coords[3] = {
  "Rstart Rend Qstart Qend Nmers Rcons Qcons Rcover Qcover Rlen Qlen Rname Qname",
  "277 928 3023 3647 198 124 125 459 454 1287 3668 pb 234R_239F",
  "277 1167 1452 601 294 185 198 647 622 1287 3000 pb 98345F_1234567R_5472F"
};

static const char* comp_details[2] = {
  "pb 234R_239F 48:-3167 49:-3166 51:-3165 52:-3164 53:-3162 54:-3161 55:-3160 57:-3157 59:-3156 60:-3155"
  " 61:-3154 62:-3153 63:-3152 64:-3151 65:-3149 66:-3148 67:-3147 70:-3145 149:-3084 [277:3023] [278:3024]"
  " [299:3047] [300:3048] [301:3049] [302:3050] [303:3051] [305:3053] [337:3082] [338:3083] [339:3084]"
  " [342:3087] [343:3088] [344:3089] [346:3091] [347:3092] [349:3094] [350:3096] [352:3098] [354:3100]"
  " [355:3101] [356:3102] [358:3103] [360:3105] [361:3106] [362:3107] [363:3108] [364:3109] [365:3110]"
  " [366:3111] [367:3112] [368:3113] [370:3115] [372:3117] [373:3118] [375:3120] [376:3121] [377:3122]"
  " [400:3143] [401:3144] [402:3145] [404:3147] [405:3148] [407:3149] [458:3192] [459:3193] [460:3194]"
  " [462:3196] [463:3197] [464:3198] [465:3199] [466:3200] [468:3202] [471:3205] [472:3206] [473:3207]"
  " [474:3208] [476:3210] [477:3211] [478:3212] [479:3213] [504:3236] [505:3237] [506:3238] [507:3239]"
  " [508:3241] [509:3242] [510:3243] [512:3245] [513:3246] [514:3247] [515:3248] [516:3249] [517:3250]"
  " [518:3251] [519:3252] [520:3253] [521:3254] [522:3255] [524:3257] [525:3258] [526:3259] [527:3260]"
  " [529:3263] [530:3264] [531:3265] [532:3266] [534:3268] [538:3271] [542:3273] [543:3274] [544:3275]"
  " [545:3276] [546:3277] [548:3279] [550:3280] [551:3281] [553:3283] [554:3284] [555:3285] [559:3289]"
  " [561:3290] [607:3338] [610:3339] [612:3341] [613:3342] [614:3343] 614:3189 [615:3344] [616:3346]"
  " [617:3347] [618:3349] [619:3350] [620:3351] [626:3355] [627:3356] [628:3357] [630:3359] [631:3360]"
  " [632:3361] [633:3362] [635:3364] [636:3365] [638:3367] [639:3368] [641:3370] [642:3371] [644:3373]"
  " [646:3375] [647:3376] [649:3378] [650:3379] [651:3380] [652:3381] [655:3384] [656:3385] [658:3387]"
  " [663:3391] [664:3392] [665:3393] [668:3396] [669:3397] [671:3399] [672:3400] [674:3402] [675:3403]"
  " [676:3404] [677:3406] [678:3407] [680:3409] [682:3410] [745:3466] [747:3468] [748:3469] [751:3473]"
  " [754:3476] [755:3477] [783:3505] [785:3506] [787:3508] [788:3509] [839:3557] [840:3560] [841:3561]"
  " [842:3562] [844:3564] [845:3565] [848:3568] [849:3569] [850:3570] [852:3572] [853:3573] [855:3575]"
  " [856:3576] [857:3577] [858:3578] [860:3580] [862:3583] [863:3584] [864:3585] [865:3586] [888:3608]"
  " [890:3610] [891:3611] [892:3612] [893:3613] [894:3614] [895:3615] [897:3617] [900:3619] [902:3621]"
  " [903:3622] [904:3623] [905:3624] [907:3626] [908:3627] [909:3628] [910:3629] [912:3631]",
  "pb 98345F_1234567R_5472F 48:1292 49:1293 51:1295 52:1296 53:1297 54:1298 55:1299 57:1301 59:1303 60:1304 61:1305 62:1306"
  " 63:1307 64:1308 65:1309 66:1310 67:1311 70:1314 149:1379 [277:-1436] [278:-1435] [299:-1413]"
  " [300:-1412] [301:-1411] [302:-1410] [303:-1409] [305:-1408] [337:-1382] [338:-1381] [339:-1379]"
  " [342:-1376] [343:-1374] [344:-1372] [346:-1370] [347:-1369] [349:-1368] [350:-1367] [352:-1365]"
  " [354:-1363] [355:-1361] [356:-1360] [358:-1359] [360:-1358] [361:-1357] [362:-1356] [363:-1353]"
  " [364:-1352] [365:-1351] [366:-1349] [367:-1348] [368:-1346] [370:-1344] [372:-1342] [373:-1340]"
  " [375:-1339] [376:-1338] [377:-1337] [400:-1318] [401:-1316] [402:-1314] [404:-1311] [405:-1310]"
  " [407:-1309] [458:-1266] [459:-1265] [460:-1264] [462:-1263] [463:-1262] [464:-1260] [465:-1259]"
  " [466:-1258] [468:-1256] [471:-1255] [472:-1254] [473:-1253] [474:-1251] [476:-1250] [477:-1249]"
  " [478:-1248] [479:-1247] [504:-1223] [505:-1221] [506:-1220] [507:-1219] [508:-1217] [509:-1216]"
  " [510:-1214] [512:-1213] [513:-1212] [514:-1211] [515:-1210] [516:-1208] [517:-1207] [518:-1206]"
  " [519:-1205] [520:-1204] [521:-1203] [522:-1202] [524:-1201] [525:-1199] [526:-1198] [527:-1197]"
  " [529:-1195] [530:-1194] [531:-1193] [532:-1192] [534:-1191] [538:-1190] [542:-1189] [543:-1188]"
  " [544:-1187] [545:-1186] [546:-1185] [548:-1183] [550:-1182] [551:-1181] [553:-1180] [554:-1177]"
  " [555:-1176] [559:-1175] [561:-1174] [607:-1123] [610:-1122] [612:-1119] [613:-1118] 614:-1272"
  " [614:-1117] [615:-1116] [616:-1114] [617:-1113] [618:-1111] [619:-1110] [620:-1109] [626:-1108]"
  " [627:-1107] [628:-1105] [630:-1103] [631:-1102] [632:-1101] [633:-1099] [635:-1098] [636:-1097]"
  " [638:-1096] [639:-1094] [641:-1093] [642:-1091] [644:-1090] [646:-1089] [647:-1085] [649:-1084]"
  " [650:-1083] [651:-1081] [652:-1080] [655:-1079] [656:-1078] [658:-1076] [663:-1075] [664:-1073]"
  " [665:-1072] [668:-1070] [669:-1069] [671:-1067] [672:-1065] [674:-1064] [675:-1062] [676:-1061]"
  " [677:-1060] [678:-1059] [680:-1056] [682:-1055] [745:-997] [747:-996] [748:-993] [751:-992] [754:-991]"
  " [755:-989] [783:-958] [785:-957] [787:-955] [788:-951] [839:-903] [840:-902] [841:-901] [842:-900]"
  " [844:-898] [845:-897] [848:-896] [849:-895] [850:-894] [852:-892] [853:-891] [855:-890] [856:-888]"
  " [857:-886] [858:-885] [860:-883] [862:-880] [863:-879] [864:-878] [865:-876] [888:-850] [890:-849]"
  " [891:-848] [892:-847] [893:-846] [894:-845] [895:-844] [897:-843] [900:-841] [902:-840] [903:-839]"
  " [904:-838] [905:-837] [907:-835] [908:-834] [909:-832] [910:-830] [912:-829] [939:-802] [941:-801]"
  " [942:-800] [946:-799] [947:-798] [949:-797] [953:-796] [954:-795] [955:-794] [956:-793] [957:-792]"
  " [958:-790] [962:-789] [963:-788] [964:-787] [965:-786] [966:-784] [967:-783] [968:-780] [970:-779]"
  " [972:-777] [974:-773] [976:-772] [977:-771] [978:-770] [980:-769] [981:-768] [1031:-716] [1034:-715]"
  " [1035:-714] [1036:-712] [1037:-711] [1038:-710] [1039:-709] [1041:-708] [1042:-707] [1043:-706]"
  " [1045:-705] [1046:-704] [1048:-702] [1049:-701] [1051:-700] [1052:-699] [1053:-696] [1054:-695]"
  " [1055:-694] [1056:-693] [1058:-692] [1059:-691] [1060:-689] [1062:-688] [1065:-687] [1066:-685]"
  " [1067:-684] [1068:-683] [1069:-682] [1094:-658] [1095:-657] [1098:-656] [1100:-655] [1101:-654]"
  " [1103:-653] [1104:-651] [1106:-650] [1107:-649] [1109:-645] [1110:-644] [1111:-642] [1113:-641]"
  " [1114:-640] [1116:-638] [1117:-637] [1118:-636] [1122:-633] [1123:-631] [1124:-630] [1125:-629]"
  " [1126:-628] [1128:-626] [1130:-625] [1131:-623] [1132:-622] [1133:-621] [1136:-618] [1138:-617]"
  " [1140:-616] [1141:-615] [1142:-614] [1143:-610] [1144:-609] [1146:-608] [1147:-607] [1148:-606]"
  " [1149:-604] [1150:-602] [1151:-601]"
};

class AlignerOutput : public ::testing::Test {
public:
  AlignerOutput() :
    sr_file("/tmp/superreads.fa"), pb_file("/tmp/pb.fa"),
    coords_file("/tmp/test.coords"), details_file("/tmp/test.details") { coords_file.do_unlink = false; details_file.do_unlink = false; }

protected:
  remove_file sr_file, pb_file, coords_file, details_file;

  void write_to(const char* path, const char* content) {
    std::ofstream os(path);
    if(!os.good())
      std::runtime_error("Failed to open fasta file");
    os.write(content, strlen(content));
  }

  virtual void SetUp() {
    write_to(sr_file.path, super_reads);
    write_to(pb_file.path, pac_bio);
  }
};

void check_file(const char* path, const char* lines[], size_t nlen, size_t header) {
  std::ifstream is(path);
  std::string line;

  for(size_t i = 0; i < header; ++i) { // All lines in header are in order
    ASSERT_TRUE(is.good());
    std::getline(is, line);
    ASSERT_STREQ(lines[i], line.c_str());
  }

  // Put remaining lines in a set to match out of order
  std::set<std::string> set_lines;
  for(size_t i = header; i < nlen; ++i)
    set_lines.insert(std::string(lines[i]));

  for(size_t i = header; i < nlen; ++i) {
    ASSERT_TRUE(is.good());
    std::getline(is, line);
    SCOPED_TRACE(::testing::Message() << "line:" << line);
    EXPECT_EQ((size_t)1, set_lines.erase(line));
  }
  EXPECT_FALSE(std::getline(is, line));
}

TEST_F(AlignerOutput, Normal) {
  mer_dna::k(17);
  mer_pos_hash_type hash(strlen(super_reads) * 2);
  frag_lists names(1);
  superread_parse(1, hash, names, sr_file.path);
  align_pb_reads(1, hash, 10, 2, 4, 10, false, false, pb_file.path, coords_file.path, details_file.path);
  check_file(coords_file.path, normal_coords, sizeof(normal_coords) / sizeof(char*), 1);
  check_file(details_file.path, normal_details, sizeof(normal_details) / sizeof(char*), 0);
} // AlignerOutput.Normal

TEST_F(AlignerOutput, Forward) {
  mer_dna::k(17);
  mer_pos_hash_type hash(strlen(super_reads) * 2);
  frag_lists names(1);
  superread_parse(1, hash, names, sr_file.path);
  align_pb_reads(1, hash, 10, 2, 4, 10, false, true, pb_file.path, coords_file.path, details_file.path);
  check_file(coords_file.path, forward_coords, sizeof(normal_coords) / sizeof(char*), 1);
  check_file(details_file.path, normal_details, sizeof(normal_details) / sizeof(char*), 0);
} // AlignerOutput.Forward


TEST_F(AlignerOutput, Compressed) {
  mer_dna::k(17);
  mer_pos_hash_type hash(strlen(super_reads) * 2);
  frag_lists names(1);
  superread_parse(1, hash, names, sr_file.path, true);
  align_pb_reads(1, hash, 10, 2, 4, 10, true, false, pb_file.path, coords_file.path, details_file.path);
  check_file(coords_file.path, comp_coords, sizeof(normal_coords) / sizeof(char*), 1);
  check_file(details_file.path, comp_details, sizeof(normal_details) / sizeof(char*), 0);

} // AlignerOutput.Compressed

} // namespace
