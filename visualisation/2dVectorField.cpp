#include <iostream>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/io/writers/VolWriter.h"
#include "DGtal/io/readers/ITKReader.h"
#include "DGtal/images/ImageSelector.h"
#include "DGtal/io/writers/ITKWriter.h"
#include "DGtal/io/boards/Board2D.h"
#include "DGtal/io/colormaps/GradientColorMap.h"

using namespace std;
using namespace DGtal;
namespace po = boost::program_options;

int main(int argc, char** argv) {

	po::options_description general_opt("Allowed options are: ");
	general_opt.add_options()
		("help,h", "display this message")
		("input,i", po::value<std::string>(), "vector field (.mha)" )
        ("fixed,f", po::value<std::string>(), "fixed image (itk format)" )
        ("arrowsize,a", po::value<int>()->default_value(1), "arrow size")
        ("spacing,s", po::value<float>()->default_value(1.0), "spacing between arrows (in pixels)")
		("output,o",  po::value<std::string>(), "output itk file" ) ;

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
	string outputFilename = vm["output"].as<std::string>();
    string fixedFilename = vm["fixed"].as<std::string>();

    int size = vm["arrowsize"].as<int>();
    float spacing = vm["spacing"].as<float>();

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
    Domain aDomain = dgtalImage.domain();

    FixedImage fixed = DGtal::ITKReader<FixedImage>::importITK(fixedFilename);

    Board2D board;
    Z2i::Point argMaxNorm = *std::max_element(aDomain.begin(), aDomain.end(), [&](const Z2i::Point& p1, const Z2i::Point& p2) {
            auto v1 = dgtalImage(p1);
            auto v2 = dgtalImage(p2);
            Z2i::RealVector v1dgtal(v1[0], v1[1]);
            Z2i::RealVector v2dgtal(v2[0], v2[1]);
            return v1dgtal.norm() < v2dgtal.norm();
        });
    float maxNorm = Z2i::RealVector(dgtalImage(argMaxNorm)[0], dgtalImage(argMaxNorm)[1]).norm();
    GradientColorMap<double> colormap( 0.0, maxNorm );
    colormap.addColor( Color( 0,   0, 0 ) );
    colormap.addColor( Color( 0,   0, 255 ) );
    colormap.addColor( Color( 0, 255, 0 ) );
    colormap.addColor( Color( 255, 0, 0 ) );
    colormap.addColor( Color( 255, 255, 0 ) );
    colormap.addColor( Color( 255, 255, 255 ) );

    std::string specificStyle =  aDomain.className() + "/Paving";
    int cpt = 0;
    std::cout << aDomain << std::endl;
	for (auto it = aDomain.begin(), ite = aDomain.end(); it != ite; ++it) {
        Z2i::Point p = *it;
         p[1] = aDomain.upperBound()[1] - p[1];
        board << CustomStyle( (*it).className(),
                              new CustomColors( Color::None,
                                                (unsigned char) fixed(*it) ) )
              << p;
    }
    for (auto it = aDomain.begin(), ite = aDomain.end(); it != ite; ++it) {
        if (cpt < spacing) {
            cpt++;
            continue;
        }
        else {
            cpt = 0;
            Z2i::Point p = *it;
            p[1] = aDomain.upperBound()[1] - p[1];
            auto vectorITK = dgtalImage(*it);
            Z2i::RealVector vector(vectorITK[0], -vectorITK[1]);
            double norm = vector.norm();
            board.setPenColor(colormap(norm));
            board.setLineWidth(0.5);
            Z2i::RealPoint destination = p + size*vector;
            board.drawArrow((float)p[0], (float)p[1], destination[0], destination[1]);
        }
	}
    board.scaleAll(0.2);
//    board.rotate(M_PI, LibBoard::Point(0,0) );
    board.saveSVG(outputFilename.c_str(), aDomain.upperBound()[0], aDomain.upperBound()[1]);
    DGtal::trace.info() << "exported" << std::endl;


	return 0;
}
