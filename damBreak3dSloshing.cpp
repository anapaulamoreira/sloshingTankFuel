/* This file is part of the Palabos library.
 *
 * Copyright (C) 2011-2017 FlowKit Sarl
 * Route d'Oron 2
 * 1010 Lausanne, Switzerland
 * E-mail contact: contact@flowkit.com
 *
 * The most recent release of Palabos can be downloaded at 
 * <http://www.palabos.org/>
 *
 * The library Palabos is free software: you can redistribute it and/or
 * modify it under the terms of the GNU Affero General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * The library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

/* The breaking dam free surface problem. This code demonstrates the basic usage of the
 * free surface module in Palabos. Surface tension and contact angles are optional. 
 */

#include "palabos3D.h"
#include "palabos3D.hh"

using namespace plb;

#define DESCRIPTOR descriptors::ForcedD3Q19Descriptor

typedef double T;


// Constante Smagorinsky para o modelo LES.
const T cSmago = 0.14;

// Dimensões físicas do sistema (em metros).
const T lx = 3.22;
const T ly = 1.0;
const T lz = 1.0;

const T rhoEmpty = T(1);
    
plint writeImagesIter   = 20;
plint getStatisticsIter = 20;

plint maxIter;
plint N;
plint nx, ny, nz;
T delta_t, delta_x;
Array<T,3> externalForce;
T nuPhys, nuLB, tau, omega, Bo, surfaceTensionLB, contactAngle;

std::string outDir;
plint obstacleCenterXYplane, obstacleLength, obstacleWidth, obstacleHeight, beginWaterReservoir, waterReservoirHeight;
plint waterLevelOne, waterLevelTwo, waterLevelThree, waterLevelFour;

// Arquivo para salvar os dados de pressão
std::ofstream out;
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*FUNCOES*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
void setupParameters() {
    delta_x = lz / N;
    nx = util::roundToInt(lx / delta_x);
    ny = util::roundToInt(ly / delta_x);
    nz = util::roundToInt(lz / delta_x);

    T gLB = 9.8 * delta_t * delta_t/delta_x; // Gravidade em unidade de lattice
    externalForce = Array<T,3>(0., 0., -gLB);
    tau            = (nuPhys*DESCRIPTOR<T>::invCs2*delta_t)/(delta_x*delta_x) + 0.5;
    omega          = 1./tau;    
    nuLB           = (tau-0.5)*DESCRIPTOR<T>::cs2; // Viscosidade em unidade de lattice
    
    surfaceTensionLB = rhoEmpty * gLB * N * N / Bo;

    obstacleCenterXYplane = util::roundToInt(0.744*N);
    obstacleLength        = util::roundToInt(0.403*N);
    obstacleWidth         = util::roundToInt(0.161*N);
    obstacleHeight        = util::roundToInt(0.161*N);
    beginWaterReservoir   = util::roundToInt((0.744+1.248)*N);
    waterReservoirHeight  = util::roundToInt(0.55*N); // 0,55 m de altura de liquido
    
    waterLevelOne   = util::roundToInt(0.496*N);
    waterLevelTwo   = util::roundToInt(2.*0.496*N);
    waterLevelThree = util::roundToInt(3.*0.496*N);
    waterLevelFour  = util::roundToInt((3.*0.496 + 1.150)*N);
}

// Especifica a condição inicial do fluido 
// (a cada célula é atribuído o sinalizador "fluido", "vazio" ou "parede").
int initialFluidFlags(plint iX, plint iY, plint iZ) {
    // Obstáculo na extremidade esquerda, que é atingido pelo fluido.
    // Se insideObstacle = false => naso havera obstaculo
    bool insideObstacle =
        iX >= obstacleCenterXYplane-obstacleWidth/2 &&
        iX <= obstacleCenterXYplane+obstacleWidth/2 &&
        iY >= ny/2-obstacleLength/2 &&
        iY <= ny/2+obstacleLength/2 &&
        iZ <= obstacleHeight+1;

    insideObstacle = false;
    if (insideObstacle) {
        return freeSurfaceFlag::wall;
    }
    else if (iX >= beginWaterReservoir && iZ <= waterReservoirHeight) {
        return freeSurfaceFlag::fluid;
    }
    else {
        return freeSurfaceFlag::empty;
    }
}

void writeResults(MultiBlockLattice3D<T,DESCRIPTOR>& lattice, MultiScalarField3D<T>& volumeFraction, plint iT)
{
    static const plint nx = lattice.getNx();
    static const plint ny = lattice.getNy();
    static const plint nz = lattice.getNz();
    Box3D slice(0, nx-1, ny/2, ny/2, 0, nz-1);

    //std::auto_ptr<MultiScalarField3D<T>> ias = computeVelocityNorm(lattice, slice);

    // Imprimir em arquivos '.gif'

    ImageWriter<T> imageWriter("leeloo");
     imageWriter.writeScaledPpm(createFileName("u", iT, 6),
                               *computeVelocityNorm(lattice, slice)); 

    imageWriter.writeScaledPpm(createFileName("rho", iT, 6),
                               *computeDensity(lattice, slice));
                   
    imageWriter.writeScaledPpm(createFileName("volumeFraction", iT, 6), *extractSubDomain(volumeFraction, slice));

    // Algoritmo para reconstruir a superfície livre e gravar um arquivo STL.
    std::vector<T> isoLevels;
    isoLevels.push_back((T) 0.5);
    typedef TriangleSet<T>::Triangle Triangle;
    std::vector<Triangle> triangles;
    isoSurfaceMarchingCube(triangles, volumeFraction, isoLevels, volumeFraction.getBoundingBox());
    TriangleSet<T>(triangles).writeBinarySTL(createFileName(outDir+"/interface", iT, 6)+".stl");

    // Imprimir em arquivo '.vti'
    VtkImageOutput3D<T> vtkOut(createFileName("Velocity", iT, 6), 1.);
    vtkOut.writeData<float>(*computeVelocityNorm(lattice), "u", 1.);
    vtkOut.writeData<3,float>(*computeVelocity(lattice), "velocity", delta_x/delta_t);
   
    // Monitorando a pressao de impacto a esquerda no fundo do tanque
    T densityData;
    float ix, iy, iz;
    float pressure;
    float fluidDensity = 1000;
    ix = 1;
    iy = ny/2;
    iz = 1;
    densityData = lattice.get(ix,iy,iz).computeDensity();
    pressure = DESCRIPTOR<T>::cs2*(densityData-1.)*delta_x*delta_x/(delta_t*delta_t)*fluidDensity;
    out << densityData << "\t\t" << pressure << std::endl;
    std::cout << "Density = " << densityData << std::endl;
    std::cout << "Pressure = " << pressure << std::endl;
}

void writeStatistics(FreeSurfaceFields3D<T,DESCRIPTOR>& fields) {
    pcout << " -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*- " << std::endl;
    T averageMass = freeSurfaceAverageMass<T,DESCRIPTOR>(fields.freeSurfaceArgs, fields.lattice.getBoundingBox());
    pcout << "Average Mass: " << averageMass  << std::endl;
    T averageDensity = freeSurfaceAverageDensity<T,DESCRIPTOR>(fields.freeSurfaceArgs, fields.lattice.getBoundingBox());
    pcout << "Average Density: " << std::setprecision(12) << averageDensity  << std::endl;

    T averageVolumeFraction = freeSurfaceAverageVolumeFraction<T,DESCRIPTOR>(fields.freeSurfaceArgs, fields.lattice.getBoundingBox());
    pcout << "Average Volume-Fraction: " << std::setprecision(12) << averageVolumeFraction  << std::endl;

    pcout << " -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*- " << std::endl;
}

//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*FUNCAO PRINCIPAL*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
int main(int argc, char **argv)
{
    out.open("pressureData.txt"); // Abre o arquivo para salvar os dados de pressão
    out << "Densidade*" << "\t" << "Pressao (Pa)" << std::endl;
    plbInit(&argc, &argv);
    global::directories().setInputDir("./");
        
    if (global::argc() != 8) {
        pcout << "Error missing some input parameter\n";
    }

    try {
        global::argv(1).read(outDir);
        global::directories().setOutputDir(outDir+"/");

        global::argv(2).read(nuPhys);
        global::argv(3).read(Bo);
        global::argv(4).read(contactAngle);
        global::argv(5).read(N);
        global::argv(6).read(delta_t);
        global::argv(7).read(maxIter);
    }
    catch(PlbIOException& except) {
        pcout << except.what() << std::endl;
        pcout << "The parameters for this program are :\n";
        pcout << "1. Output directory name.\n";
        pcout << "2. kinematic viscosity in physical Units (m^2/s) .\n";
        pcout << "3. Bond number (Bo = rho * g * L^2 / gamma).\n";
        pcout << "4. Contact angle (in degrees).\n";
        pcout << "5. number of lattice nodes for lz .\n";
        pcout << "6. delta_t .\n"; 
        pcout << "7. maxIter .\n";
        pcout << "Reasonable parameters on a desktop computer are: " << (std::string)global::argv(0) << " tmp 1.e-5 100 80.0 40 1.e-3 80000\n";
        pcout << "Reasonable parameters on a parallel machine are: " << (std::string)global::argv(0) << " tmp 1.e-6 100 80.0 100 1.e-4 80000\n";
        exit (EXIT_FAILURE);
    }
    
    setupParameters();
    
    pcout << "delta_t= " << delta_t << std::endl;
    pcout << "delta_x= " << delta_x << std::endl;
    pcout << "delta_t*delta_t/delta_x= " << delta_t*delta_t/delta_x << std::endl;
    pcout << "externalForce= " << externalForce[2] << std::endl;
    pcout << "relaxation time= " << tau << std::endl;
    pcout << "omega= " << omega << std::endl;
    pcout << "kinematic viscosity physical units = " << nuPhys << std::endl;
    pcout << "kinematic viscosity lattice units= " << nuLB << std::endl;
    
    global::timer("initialization").start();
    
    // Um MultiBlockDistribution3D descreve o arranjo de blocos atômicos dentro de um multi-bloco. 
    // Cada um deles obtém um ID inteiro, esse ID é usada em programas paralelos para associar um bloco a um processo MPI.
    // Para nao haver desperdício de memória, pois ha regioes que não são utilizadas, ao inves de representar esse
    // problema por uma estrutura de matriz regular, eh utilizado o domínio esparso.
    SparseBlockStructure3D blockStructure(createRegularDistribution3D(nx, ny, nz));

    Dynamics<T,DESCRIPTOR>* dynamics
        = new SmagorinskyBGKdynamics<T,DESCRIPTOR>(omega, cSmago);

    // Se surfaceTensionLB for 0, o algoritmo de tensão superficial será desativado.
    // Se contactAngle for menor que 0, o algoritmo de ângulo de contato será desativado.
    FreeSurfaceFields3D<T,DESCRIPTOR> fields( blockStructure, dynamics->clone(), rhoEmpty,
                                              surfaceTensionLB, contactAngle, externalForce );

    // A função integrateProcessingFunctional é usada para adicionar um processador de dados a um bloco e 
    // executá-lo periodicamente a cada ciclo de iteração. Aqui por conta das condições de contorno não locais.
    //integrateProcessingFunctional(new ShortenBounceBack3D<T,DESCRIPTOR>, fields.lattice.getBoundingBox(), fields.freeSurfaceArgs, 0);

    // Define todas as células da parede externa para "parede" (aqui, as células em massa também são definidas 
    // para "parede", mas não importa, porque eles são sobrescritos na próxima linha).
    setToConstant(fields.flag, fields.flag.getBoundingBox(), (int)freeSurfaceFlag::wall);

    // No volume (todos, exceto na camada externa), inicializa os sinalizadores conforme especificado pela
    // funcao "initialFluidFlags".
    setToFunction(fields.flag, fields.flag.getBoundingBox().enlarge(-1), initialFluidFlags);
    
    fields.defaultInitialize();

    pcout << "Time spent for setting up lattices: "
          << global::timer("initialization").stop() << std::endl;
    T lastIterationTime = T();

//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*LACO PRINCIPAL*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
    for (plint iT = 0; iT <= maxIter; ++iT) {
        global::timer("iteration").restart();
        
        T sum_of_mass_matrix = T();
        T lost_mass = T();
        if (iT % getStatisticsIter==0) {
            pcout << std::endl;
            pcout << "ITERATION = " << iT << std::endl;
            pcout << "Time of last iteration is " << lastIterationTime << " seconds" << std::endl;
            writeStatistics(fields);
            sum_of_mass_matrix = fields.lattice.getInternalStatistics().getSum(0);
            pcout << "Sum of mass matrix: " << sum_of_mass_matrix << std::endl;
            lost_mass = fields.lattice.getInternalStatistics().getSum(1);
            pcout << "Lost mass: " << lost_mass << std::endl;
            pcout << "Total mass: " << sum_of_mass_matrix + lost_mass << std::endl;
            pcout << "Interface cells: " << fields.lattice.getInternalStatistics().getIntSum(0) << std::endl;
        }

        if (iT % writeImagesIter == 0) {
            global::timer("images").start();
            writeResults(fields.lattice, fields.volumeFraction, iT);
            pcout << "Total time spent for writing images: "
                << global::timer("images").stop() << std::endl;
        }                           

        // Isso inclui o ciclo de colisao e adveccao e todas as operacoes de superficie livre. 
        // Mais detalhes na observacao abaixo
        fields.lattice.executeInternalProcessors();
        fields.lattice.evaluateStatistics();
        fields.lattice.incrementTime();

        lastIterationTime = global::timer("iteration").stop();
    }
    out.close(); // Fechamento do arquivo
}