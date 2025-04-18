

// core headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/pose/extra_pose_info_util.hh>
#include <devel/init.hh>
#include <basic/options/option_macros.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
//#include <protocols/relax_protocols.hh>
#include <basic/Tracer.hh>
#include <core/scoring/hbonds/hbonds.hh>
#include <core/scoring/hbonds/HBondSet.hh>


// utility headers
#include <utility/vector1.hh>
#include <utility/file/FileName.fwd.hh>

// c++ headers
#include <string>

#include <core/scoring/hbonds/HBondDatabase.hh>
#include <core/scoring/hbonds/HBondOptions.hh>

#include <protocols/jd2/JobDistributor.hh>

#include <protocols/moves/Mover.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/moves/MoverContainer.fwd.hh>
#include <protocols/symmetry/SetupForSymmetryMover.hh>


// namespaces
using namespace core;
using namespace conformation;
using namespace chemical;
using namespace scoring;
using namespace hbonds;
using namespace pose;
using namespace ObjexxFCL;
using namespace basic;
using namespace basic::options;
using namespace basic::options::OptionKeys;
using basic::Error;
using basic::Warning;
using utility::file::FileName;

OPT_1GRP_KEY(Boolean, min, scoreonly)


static basic::Tracer TR("apps.xtaldocking.hbonds");

class CrystHbondsReporter : public protocols::moves::Mover {
public:
	CrystHbondsReporter()= default;

	void apply( core::pose::Pose & pose) override {

		HBondDatabaseCOP hb_database( HBondDatabase::get_database());
		HBondOptions hboptions;
		HBondSet set1;
		utility::vector1< HBondCOP > hbonds;

		//pose.dump_pdb("cryst_test_counthbonds.pdb");

		ScoreFunctionOP scfxn( get_score_function() );
		pose.update_residue_neighbors();
		fill_hbond_set( pose, false, set1, false );
		//note that the list of HBonds is indexed starting at 1
		hbonds = set1.residue_hbonds(1);
		core::Size nhbonds( hbonds.size() );
		core::pose::setPoseExtraScore( pose, "nhbonds", nhbonds);
		(*scfxn)(pose);

	}

	std::string get_name() const override {
		return "CrystHbondsReporter";
	}
};



int
main( int argc, char * argv [] ) {
	using namespace protocols::moves;
	using namespace protocols;
	using namespace protocols::jd2;
	using namespace protocols::symmetry;

	try {
		NEW_OPT(min::scoreonly, "scoresonly?", false);

		devel::init(argc, argv);

		SequenceMoverOP seq( new SequenceMover() );

		seq->add_mover( utility::pointer::make_shared< SetupForSymmetryMover >() );
		seq->add_mover( utility::pointer::make_shared< CrystHbondsReporter >() );

		protocols::jd2::JobDistributor::get_instance()->go( seq );
	} catch (utility::excn::Exception const & e ) {
		e.display();
		return -1;
	}
	return 0;
}
