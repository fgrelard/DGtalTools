/// @file dec/exampleDiscreteExteriorCalculusSolve.cpp
/**
   Example of primal and dual Helmoltz decomposition in 2D and 3D using Discrete Exterior Calculus.
   @see \ref sectDECHelmoltzProblem
   \image  html  solve_2d_primal_decomposition_calculusSmall.png "Primal Helmoltz decomposition harmonic component."
   \example dec/exampleDiscreteExteriorCalculusSolve.cpp
**/

#include <string>

#include <QApplication>

#include "DECExamplesCommon.h"

#include <iostream>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

// always include EigenSupport.h before any other Eigen headers
#include "DGtal/math/linalg/EigenSupport.h"
#include "DGtal/dec/DiscreteExteriorCalculus.h"
#include "DGtal/dec/DiscreteExteriorCalculusSolver.h"
#include "DGtal/dec/DiscreteExteriorCalculusFactory.h"
#include "DGtal/base/Common.h"
#include "DGtal/io/viewers/Viewer3D.h"
#include "DGtal/io/boards/Board2D.h"
#include "DGtal/io/readers/GenericReader.h"

using namespace DGtal;
using namespace std;
namespace po = boost::program_options;



void solve2d_vectorfield_decomposition(std::string inputFilename) {
    trace.beginBlock("2d discrete exterior calculus solve primal helmoltz decomposition");
    typedef Z2i::Space Space;
	typedef Z2i::KSpace KSpace;
	typedef HyperRectDomain<Space> Domain;
    typedef itk::Vector<double, 2> PixelType;
    typedef itk::Image<PixelType, 2> ImageITK;
    typedef ImageContainerByITKImage<Domain, float> FixedImage;
    typedef ImageContainerByITKImage<Domain, PixelType> ImageBridge;
    typedef itk::ImageFileReader<ImageITK> ITKReader;

    ITKReader::Pointer reader = ITKReader::New();
    reader->SetFileName(inputFilename);
    reader->Update();
    ImageITK::Pointer image = reader->GetOutput();

    ImageBridge dgtalImage(image);
    Z2i::Domain aDomain = dgtalImage.domain();


    // create discrete exterior calculus from set
    typedef DiscreteExteriorCalculus<2, 2, EigenLinearAlgebraBackend> Calculus;
    // Calculus calculus;

    typedef DiscreteExteriorCalculusFactory<EigenLinearAlgebraBackend> CalculusFactory;
    Z2i::DigitalSet aSet(aDomain);
    aSet.insert(aDomain.begin(), aDomain.end());
    Calculus calculus = CalculusFactory::createFromDigitalSet(aSet);
    // Calculus calculus = CalculusFactory::createFromDigitalSet(generateDoubleRingSet(domain));
    // trace.info() << calculus << endl;
    // choose linear solver
    typedef EigenLinearAlgebraBackend::SolverSparseQR LinearAlgebraSolver;

    //! [2d_dual_decomposition_operator_definition]
    const Calculus::DualDerivative0 d0 = calculus.derivative<0, DUAL>();
    const Calculus::DualDerivative1 d1 = calculus.derivative<1, DUAL>();
    const Calculus::PrimalDerivative0 d0p = calculus.derivative<0, PRIMAL>();
    const Calculus::PrimalDerivative1 d1p = calculus.derivative<1, PRIMAL>();
    const Calculus::DualHodge1 h1 = calculus.hodge<1, DUAL>();
    const Calculus::DualHodge2 h2 = calculus.hodge<2, DUAL>();
    const Calculus::PrimalHodge1 h1p = calculus.hodge<1, PRIMAL>();
    const Calculus::PrimalHodge2 h2p = calculus.hodge<2, PRIMAL>();
    const LinearOperator<Calculus, 1, DUAL, 0, DUAL> ad1 = h2p * d1p * h1;
    const LinearOperator<Calculus, 2, DUAL, 1, DUAL> ad2 = h1p * d0p * h2;
    //! [2d_dual_decomposition_operator_definition]

    //! [2d_dual_decomposition_input_field_definition]
    Calculus::DualVectorField input_vector_field(calculus);
    int ii = 0;
    for (auto it = aDomain.begin(), ite = aDomain.end(); it != ite; ++it) {
        Z2i::Point p = *it;
        p[1] = aDomain.upperBound()[1] - p[1];
        auto vectorITK = dgtalImage(*it);
        Z2i::RealVector vector(vectorITK[0], -vectorITK[1]);
        Z2i::RealVector vectorNorm = vector.getNormalized();
        input_vector_field.setVector(ii, vector);
        ii++;
    }
    // for (Calculus::Index ii=0; ii<input_vector_field.length(); ii++)
    // {
    //     const Z2i::RealPoint cell_center = Z2i::RealPoint(input_vector_field.getSCell(ii).preCell().coordinates)/2.;
    //     input_vector_field.myCoordinates(ii, 0) = cos(-.5*cell_center[0]+ .3*cell_center[1]);
    //     input_vector_field.myCoordinates(ii, 1) = cos(.4*cell_center[0]+ .8*cell_center[1]);
    //     trace.info() << input_vector_field.myCoordinates(ii, 0) << std::endl;
    // }
    trace.info() << input_vector_field << endl;

    const Calculus::DualForm1 input_one_form = calculus.flat(input_vector_field);
    const Calculus::DualForm0 input_one_form_anti_derivated = ad1 * input_one_form;
    const Calculus::DualForm2 input_one_form_derivated = d1 * input_one_form;
    //! [2d_dual_decomposition_input_field_definition]

    {
        Board2D board;
        board.setLineWidth(0.1);

        // board << aDomain;
        // board << calculus;
        board << CustomStyle("KForm", new KFormStyle2D(-1, 1));
        // board << input_one_form;
        board << CustomStyle("VectorField", new VectorFieldStyle2D(.75));
        int i = 0;
        for (auto it = aDomain.begin(), ite = aDomain.end(); it != ite; ++it) {
            auto vector = input_vector_field.getVector(i);
            Z2i::Point p = *it;
            Z2i::RealPoint destination = p + vector;
            board.drawArrow((float)p[0], (float)p[1], destination[0], destination[1]);
            i++;
        }
        // board << input_vector_field;

        board.scaleAll(0.1);
        board.saveSVG("calculus.svg");
    }


    Calculus::DualForm0 solution_curl_free(calculus);
    { // solve curl free problem
        trace.beginBlock("solving curl free component");

        //! [2d_dual_decomposition_curl_free_solve]
        typedef DiscreteExteriorCalculusSolver<Calculus, LinearAlgebraSolver, 0, DUAL, 0, DUAL> Solver;
        Solver solver;
        solver.compute(ad1 * d0);
        solution_curl_free = solver.solve(input_one_form_anti_derivated);
        //! [2d_dual_decomposition_curl_free_solve]

        trace.info() << solver.isValid() << " " << solver.myLinearAlgebraSolver.info() << endl;
        trace.info() << "min=" << solution_curl_free.myContainer.minCoeff() << " max=" << solution_curl_free.myContainer.maxCoeff() << endl;
        trace.endBlock();
    }

    {
        Board2D board;
        board << aDomain;
        board << calculus;
        // board << solution_curl_free;
        board << CustomStyle("VectorField", new VectorFieldStyle2D(.75));
        board << calculus.sharp(d0*solution_curl_free);
        board.saveSVG("curl_free.svg");
    }

    Calculus::DualForm2 solution_div_free(calculus);
    { // solve divergence free problem
        trace.beginBlock("solving divergence free component");

        //! [2d_dual_decomposition_div_free_solve]
        typedef DiscreteExteriorCalculusSolver<Calculus, LinearAlgebraSolver, 2, DUAL, 2, DUAL> Solver;
        Solver solver;
        solver.compute(d1 * ad2);
        solution_div_free = solver.solve(input_one_form_derivated);
        //! [2d_dual_decomposition_div_free_solve]

        trace.info() << solver.isValid() << " " << solver.myLinearAlgebraSolver.info() << endl;
        trace.info() << "min=" << solution_div_free.myContainer.minCoeff() << " max=" << solution_div_free.myContainer.maxCoeff() << endl;
        trace.endBlock();
    }

    {
        Board2D board;
        board << aDomain;
        board << calculus;
        // board << solution_div_free;
        board << CustomStyle("VectorField", new VectorFieldStyle2D(1.5));
        board << calculus.sharp(ad2*solution_div_free);
        board.saveSVG("div_free.svg");
    }

    //! [2d_dual_decomposition_solution]
    const Calculus::DualForm1 solution_harmonic = input_one_form - d0*solution_curl_free - ad2*solution_div_free;
    //! [2d_dual_decomposition_solution]
    trace.info() << "min=" << solution_harmonic.myContainer.minCoeff() << " max=" << solution_harmonic.myContainer.maxCoeff() << endl;

    {
        Board2D board;
        board << aDomain;
        board << calculus;
        // board << solution_harmonic;
        board << CustomStyle("VectorField", new VectorFieldStyle2D(20));
        board << calculus.sharp(solution_harmonic);
        board.saveSVG("harmonic.svg");
    }

    trace.endBlock();
}


int main(int argc, char* argv[])
{
    po::options_description general_opt("Allowed options are: ");
	general_opt.add_options()
		("help,h", "display this message")
		("input,i", po::value<std::string>(), "vector field (.mha)" )
        // ("fixed,f", po::value<std::string>(), "fixed image (itk format)" )
        // ("arrowsize,a", po::value<int>()->default_value(1), "arrow size")
        // ("spacing,s", po::value<float>()->default_value(1.0), "spacing between arrows (in pixels)")
		// ("output,o",  po::value<std::string>(), "output itk file" )
        ;

	bool parseOK=true;
	po::variables_map vm;
	try{
		po::store(po::parse_command_line(argc, argv, general_opt), vm);
	}catch(const std::exception& ex){
		parseOK=false;
		trace.info()<< "Error checking program options: "<< ex.what()<< endl;
	}
	po::notify(vm);
	if( !parseOK || vm.count("help")||argc<=1)
	{
		std::cout << "Usage: " << argv[0] << " [input]\n"
				  << "Display volume file as a voxel set by using QGLviewer"<< endl
				  << general_opt << "\n";
		return 0;
	}
	if(!vm.count("input"))
	{
		trace.error() << " The file name was not defined" << endl;
		return 0;
	}
	string inputFilename = vm["input"].as<std::string>();
	// string outputFilename = vm["output"].as<std::string>();
    // string fixedFilename = vm["fixed"].as<std::string>();

    // int size = vm["arrowsize"].as<int>();
    // float spacing = vm["spacing"].as<float>();

    // QApplication app(argc,argv);
    solve2d_vectorfield_decomposition(inputFilename);
    return 0;
    // solve2d_laplace();
    // solve2d_dual_decomposition();
    // solve2d_primal_decomposition();
    // solve3d_decomposition();

    // return app.exec();
}
