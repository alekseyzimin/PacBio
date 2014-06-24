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

static const char* normal_coords[4] = {
  "Rstart Rend Qstart Qend Nmers Rcons Qcons Rcover Qcover Rlen Qlen Stretch Offset Err Rname Qname",
  "303 877 3051 3597 90 79 79 263 262 1287 3668 1.05042 -2898.38 1.72222 pb 1R_3F",
  "303 1150 1420 613 107 93 93 328 327 1287 3000 -1.04824 1794.59 1.4486 pb 5F_4R_2F",
  "303 1150 1280 473 107 93 93 328 327 1287 2800 -1.04824 1647.84 1.4486 pb 7R_2F"
};

static const char* forward_coords[4] = {
  "Rstart Rend Qstart Qend Nmers Rcons Qcons Rcover Qcover Rlen Qlen Stretch Offset Err Rname Qname",
  "303 877 3051 3597 90 79 79 263 262 1287 3668 1.05042 -2898.38 1.72222 pb 1R_3F",
  "303 1150 1581 2388 107 93 93 328 327 1287 3000 1.04824 -1350.17 1.4486 pb 2R_4F_5R",
  "303 1150 1521 2328 107 93 93 328 327 1287 2800 1.04824 -1287.28 1.4486 pb 2R_7F"
};

static const char* merinfo_coords[4] = {
  "Rstart Rend Qstart Qend Nmers Rcons Qcons Rcover Qcover Rlen Qlen Stretch Offset Err Rname Qname",
  "303 877 3051 3597 90 79 79 263 262 1287 3668 1.05042 -2898.38 1.72222 pb 1R_3F 0:0 0:0 90:262",
  "303 1150 1581 2388 107 93 93 328 327 1287 3000 1.04824 -1350.17 1.4486 pb 2R_4F_5R 0:0 0:0 82:222 5:21 30:126",
  "303 1150 1521 2328 107 93 93 328 327 1287 2800 1.04824 -1287.28 1.4486 pb 2R_7F 0:0 0:0 107:327"
};

static const char* normal_details[3] = {
  "pb 1R_3F 56:-3155 57:-3154 58:-3153 59:-3152 60:-3151 61:-3150 62:-3149 63:-3148 64:-3147 65:-3146"
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
  "pb 5F_4R_2F 56:1300 57:1301 58:1302 59:1303 60:1304 61:1305 62:1306 63:1307 64:1308 65:1309 66:1310 67:1311"
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
  " [1133:-614] [1134:-613]",
  "pb 7R_2F 56:1160 57:1161 58:1162 59:1163 60:1164 61:1165 62:1166 63:1167 64:1168 65:1169 66:1170 67:1171"
  " 68:1172 69:1173 70:1174 71:1175 [303:-1264] [304:-1263] [305:-1262] [306:-1261] [307:-1260] [345:-1225]"
  " [346:-1224] [347:-1223] [348:-1222] [349:-1221] [373:-1197] [374:-1196] [375:-1195] [376:-1194]"
  " [377:-1193] [403:-1169] [404:-1168] [405:-1167] [455:-1126] [456:-1125] [457:-1124] [458:-1123]"
  " [459:-1122] [460:-1121] [461:-1120] [462:-1119] [463:-1118] [464:-1117] [465:-1116] [466:-1115]"
  " [467:-1114] [468:-1113] [469:-1112] [470:-1111] [471:-1110] [472:-1109] [473:-1108] [474:-1107]"
  " [475:-1106] [476:-1105] [477:-1104] [478:-1103] [479:-1102] [480:-1101] [481:-1100] [523:-1059]"
  " [524:-1058] [525:-1057] [526:-1056] [527:-1055] [528:-1054] [638:-948] [639:-947] [640:-946]"
  " [641:-945] [642:-944] [643:-943] [644:-942] [645:-941] [646:-940] [647:-939] [648:-938]"
  " [649:-937] [650:-936] [651:-935] [652:-934] [653:-933] [654:-932] [655:-931] [656:-930]"
  " [657:-929] [658:-928] [659:-927] [660:-926] [661:-925] [675:-912] [676:-911] [746:-848]"
  " [747:-847] [748:-846] [749:-845] [750:-844] [813:-782] [855:-740] [856:-739] [857:-738] [858:-737]"
  " [859:-736] [860:-735] [861:-734] [959:-642] [960:-641] [961:-640] [962:-639] [963:-638] [964:-637]"
  " [965:-636] [966:-635] [967:-634] [1043:-562] [1044:-561] [1045:-560] [1046:-559] [1131:-476] [1132:-475]"
  " [1133:-474] [1134:-473]"
};

// Unitig lengths (0 and 6 are not used)
static const int unitig_lengths[8] = {
  0, 3000, 1043, 733, 1043, 1044, 0, 1822
};

class AlignerOutput : public ::testing::Test {
public:
  AlignerOutput() :
    coords_file(".aligner_test.coords"), details_file(".aligner_test.details") { coords_file.do_unlink = false; details_file.do_unlink = false; }

  static constexpr const char* pb_file = "test_pacbio.fa";
  static constexpr const char* sr_file = "test_super_reads.fa";

protected:
  virtual void SetUp() {
    std::ifstream file(sr_file);
    file.seekg(0, std::ios::end);
    super_read_approx_len = file.tellg();
  }

  size_t super_read_approx_len;
  remove_file coords_file, details_file;
};

void check_file(const char* path, const char* lines[], size_t nlen, size_t header) {
  std::ifstream is(path);
  std::string line;

  SCOPED_TRACE(::testing::Message() << "Check file " << path);
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
  mer_pos_hash_type hash(super_read_approx_len * 2);
  frag_lists names(1);
  superread_parse(1, hash, names, sr_file);
  align_pb_reads(1, hash, 10, 2, false, pb_file, coords_file.path, details_file.path);
  check_file(coords_file.path, normal_coords, sizeof(normal_coords) / sizeof(char*), 1);
  check_file(details_file.path, normal_details, sizeof(normal_details) / sizeof(char*), 0);
} // AlignerOutput.Normal

TEST_F(AlignerOutput, Forward) {
  mer_dna::k(17);
  mer_pos_hash_type hash(super_read_approx_len * 2);
  frag_lists names(1);
  superread_parse(1, hash, names, sr_file);
  align_pb_reads(1, hash, 10, 2, true, pb_file, coords_file.path, details_file.path);
  check_file(coords_file.path, forward_coords, sizeof(forward_coords) / sizeof(char*), 1);
  check_file(details_file.path, normal_details, sizeof(normal_details) / sizeof(char*), 0);
} // AlignerOutput.Forward


TEST_F(AlignerOutput, MerInfo) {
  mer_dna::k(17);
  const int mer_len = 65; // mer len for k-unitigs
  mer_pos_hash_type hash(super_read_approx_len * 2);
  frag_lists names(1);
  superread_parse(1, hash, names, sr_file);
  align_pb_reads(1, hash, 10, 2, true, pb_file, coords_file.path, details_file.path,
                 unitig_lengths, sizeof(unitig_lengths) / sizeof(int), mer_len);
  check_file(coords_file.path, merinfo_coords, sizeof(merinfo_coords) / sizeof(char*), 1);
  check_file(details_file.path, normal_details, sizeof(normal_details) / sizeof(char*), 0);
} // AlignerOutput.Duplicated

TEST_F(AlignerOutput, FilterMer) {
  mer_dna::k(17);
  const int mer_len = 65; // mer len for k-unitigs
  mer_pos_hash_type hash(super_read_approx_len * 2);
  frag_lists names(1);
  superread_parse(1, hash, names, sr_file);
  align_pb_reads(1, hash, 10, 2, true, pb_file, coords_file.path, details_file.path,
                 unitig_lengths, sizeof(unitig_lengths) / sizeof(int), mer_len,
                 0.09);
  check_file(coords_file.path, merinfo_coords, 2, 1);
} // AlignerOutput.Filter

TEST_F(AlignerOutput, FilterBases) {
  mer_dna::k(17);
  const int mer_len = 65; // mer len for k-unitigs
  mer_pos_hash_type hash(super_read_approx_len * 2);
  frag_lists names(1);
  superread_parse(1, hash, names, sr_file);
  align_pb_reads(1, hash, 10, 2, true, pb_file, coords_file.path, details_file.path,
                 unitig_lengths, sizeof(unitig_lengths) / sizeof(int), mer_len,
                 0.0, 0.27);
  check_file(coords_file.path, merinfo_coords, 2, 1);
} // AlignerOutput.FilterBases

} // namespace
