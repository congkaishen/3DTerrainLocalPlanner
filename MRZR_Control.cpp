



// =============================================================================
// PROJECT CHRONO - http://projectchrono.org
//
// Copyright (c) 2014 projectchrono.org
// All rights reserved.
//
// Use of this source code is governed by a BSD-style license that can be found
// in the LICENSE file at the top level of the distribution and at
// http://projectchrono.org/license-chrono.txt.
//
// =============================================================================
// Authors: Radu Serban
// =============================================================================
//
// Main driver function for a vehicle specified through JSON files.
//
// The vehicle reference frame has Z up, X towards the front of the vehicle, and
// Y pointing to the left.
//
// =============================================================================
#include "chrono/physics/ChBodyEasy.h"
#include "chrono/solver/ChIterativeSolverLS.h"

#include "chrono_vehicle/ChConfigVehicle.h"
#include "chrono_vehicle/ChVehicleModelData.h"
#include "chrono_vehicle/terrain/RigidTerrain.h"
#include "chrono_vehicle/utils/ChUtilsJSON.h"
#include "chrono_models/vehicle/hmmwv/HMMWV.h"
#include <chrono>
#include <thread>
#include "chrono_vehicle/wheeled_vehicle/vehicle/WheeledVehicle.h"
#include "chrono_vehicle/wheeled_vehicle/vehicle/WheeledTrailer.h"
#include "chrono_vehicle/utils/ChSpeedController.h"

#include "chrono_thirdparty/filesystem/path.h"

#ifdef CHRONO_IRRLICHT
#include "chrono_vehicle/driver/ChInteractiveDriverIRR.h"
#include "chrono_vehicle/wheeled_vehicle/ChWheeledVehicleVisualSystemIrrlicht.h"
using namespace chrono::irrlicht;
#endif
#include <chrono>
#include <thread>
#include <string>

using namespace chrono;
using namespace chrono::vehicle;
using namespace chrono::vehicle::hmmwv;



// =============================================================================
// Specification of a vehicle model from JSON files
// Available models:
//    MRZR

class Vehicle_Model {
public:
    virtual std::string ModelName() const = 0;
    virtual std::string VehicleJSON() const = 0;
    virtual std::string TireJSON() const = 0;
    virtual std::string PowertrainJSON() const = 0;
    virtual double CameraDistance() const = 0;
    virtual ChContactMethod ContactMethod() const = 0;
};

class Polaris_Model : public Vehicle_Model {
public:
    virtual std::string ModelName() const override { return "Polaris"; }
    virtual std::string VehicleJSON() const override {
        return "Polaris/Polaris.json";
    }
    virtual std::string TireJSON() const override {
        //return "mrzr/JSON_new/tire/MRZR_Pac02Tire.json";
        return "Polaris/Polaris_Pac02Tire.json";
    }
    virtual std::string PowertrainJSON() const override {
        return "Polaris/Polaris_SimplePowertrain.json";
    }
    virtual double CameraDistance() const override { return 6.0; }
    virtual ChContactMethod ContactMethod() const { return ChContactMethod::SMC; }
};

class MRZR_Model : public Vehicle_Model {
public:
    virtual std::string ModelName() const override { return "MRZR"; }
    virtual std::string VehicleJSON() const override {
        return "mrzr/JSON_new/vehicle/MRZR.json";
    }
    virtual std::string TireJSON() const override {
        return "mrzr/JSON_new/tire/MRZR_Pac02Tire.json";
    }
    virtual std::string PowertrainJSON() const override {
        return "mrzr/JSON_new/powertrain/MRZR_SimplePowertrain.json";
    }
    virtual double CameraDistance() const override { return 6.0; }
    virtual ChContactMethod ContactMethod() const { return ChContactMethod::SMC; }
};

class HMMWV_Model : public Vehicle_Model {
public:
    virtual std::string ModelName() const override { return "HMMWV"; }
    virtual std::string VehicleJSON() const override {
        return "hmmwv/vehicle/HMMWV_Vehicle_mapShock.json";
        ////return "hmmwv/vehicle/HMMWV_Vehicle.json";
        ////return "hmmwv/vehicle/HMMWV_Vehicle_bushings.json";
        ////return "hmmwv/vehicle/HMMWV_Vehicle_4WD.json";
    }
    virtual std::string TireJSON() const override {
        ////return "hmmwv/tire/HMMWV_RigidTire.json";
        ////return "hmmwv/tire/HMMWV_FialaTire.json";
        //return "hmmwv/tire/HMMWV_TMeasyTire.json";
        return "mrzr/JSON_new/tire/MRZR_Pac02Tire.json";
        ////return "hmmwv/tire/HMMWV_Pac89Tire.json";
        ////return "hmmwv/tire/HMMWV_Pac02Tire.json";
    }
    virtual std::string PowertrainJSON() const override {
        return "hmmwv/powertrain/HMMWV_ShaftsPowertrain.json";
        ////return "hmmwv/powertrain/HMMWV_SimpleCVTPowertrain.json";
        ////return "hmmwv/powertrain/HMMWV_SimplePowertrain.json";
    }
    virtual double CameraDistance() const override { return 6.0; }
    virtual ChContactMethod ContactMethod() const { return ChContactMethod::SMC; }
};

// =============================================================================


// Run-time visualization system (IRRLICHT or VSG)
ChVisualSystem::Type vis_type = ChVisualSystem::Type::IRRLICHT;

// Current vehicle model selection
auto vehicle_model = MRZR_Model();
//auto vehicle_model = Polaris_Model();
//auto vehicle_model = HMMWV_Model();

// JSON files for terrain
std::string rigidterrain_file("terrain/RigidPlane.json");
////std::string rigidterrain_file("terrain/RigidMesh.json");
////std::string rigidterrain_file("terrain/RigidHeightMap.json");
//std::string rigidterrain_file("terrain/RigidSlope10.json");
////std::string rigidterrain_file("terrain/RigidSlope20.json");

// Simulation step size
double step_size = 1e-3;
// Output directory
const std::string out_dir = GetChronoOutputPath() + "WHEELED_JSON";
bool img_output = false;
const std::string img_dir = GetChronoOutputPath() + "/MRZR";
double render_step_size = 1.0 / 60;

double WAIT_TIME = 1.0;
double initial_state[4] = { 0.0, 0.0,  0.2 , 0};
double goal[2] = { 110.0, 0.0 };
double friction = 0.8;
double slope = 0.3;

//std::string filenamer(double friction, double YAW_ANG, double speed_cmd, double sa_mag, double sa_freq) {
//    std::string num_text = std::to_string(friction);
//    std::string rounded = num_text.substr(0, num_text.find(".") + 2);
//    std::string FILENAME = rounded + "_";
//
//    num_text = std::to_string(YAW_ANG);
//    rounded = num_text.substr(0, num_text.find(".") + 3);
//    FILENAME = FILENAME + rounded + "_";
//
//    num_text = std::to_string(speed_cmd);
//    rounded = num_text.substr(0, num_text.find(".") + 2);
//    FILENAME = FILENAME + rounded + "_";
//
//    num_text = std::to_string(sa_mag);
//    rounded = num_text.substr(0, num_text.find(".") + 3);
//    FILENAME = FILENAME + rounded + "_";
//
//    num_text = std::to_string(sa_freq);
//    rounded = num_text.substr(0, num_text.find(".") + 2);
//    FILENAME = FILENAME + rounded + "_";
//
//    return FILENAME;
//}

void run_sim() {
    GetLog() << "Copyright (c) 2017 projectchrono.org\nChrono version: " << CHRONO_VERSION << "\n\n";
    SetChronoDataPath(CHRONO_DATA_DIR);
    chrono::vehicle::SetDataPath("E:/workspace/chrono_8/data/vehicle/");

    // Create the vehicle system
    WheeledVehicle vehicle(vehicle::GetDataFile(vehicle_model.VehicleJSON()), vehicle_model.ContactMethod());
    vehicle.Initialize(ChCoordsys<>(ChVector<>(initial_state[0], initial_state[1], initial_state[2]), Q_from_AngX(slope)));
    vehicle.GetChassis()->SetFixed(false);
    vehicle.SetChassisVisualizationType(VisualizationType::MESH);
    vehicle.SetChassisRearVisualizationType(VisualizationType::MESH);
    vehicle.SetSuspensionVisualizationType(VisualizationType::MESH);
    vehicle.SetSteeringVisualizationType(VisualizationType::MESH);
    vehicle.SetWheelVisualizationType(VisualizationType::MESH);
    auto powertrain = ReadPowertrainJSON(vehicle::GetDataFile(vehicle_model.PowertrainJSON()));
    vehicle.InitializePowertrain(powertrain);
    vehicle.LockCentralDifferential(0, true);
    vehicle.LockCentralDifferential(1, true);
    for (auto& axle : vehicle.GetAxles()) {
        for (auto& wheel : axle->GetWheels()) {
            auto tire = ReadTireJSON(vehicle::GetDataFile(vehicle_model.TireJSON()));
            vehicle.InitializeTire(tire, wheel, VisualizationType::MESH);
        }
    }
    //// Containing system
    //auto system = vehicle.GetSystem();
    //RigidTerrain terrain(system);
    //auto patch1_mat = chrono_types::make_shared<ChMaterialSurfaceSMC>();
    //patch1_mat->SetFriction(friction);
    //patch1_mat->SetRestitution(0.01f);
    //auto patch1 = terrain.AddPatch(patch1_mat, ChCoordsys<>(ChVector<>(0, 0, 0), Q_from_AngX(0.0)),
    //    vehicle::GetDataFile("terrain/meshes/challenge_9.obj"));
    //patch1->SetColor(ChColor(1.0f, 1.0f, 1.0f));
    //patch1->SetTexture(vehicle::GetDataFile("terrain/textures/tile4.jpg"), 200, 200);
    //terrain.Initialize();


    auto system = vehicle.GetSystem();

    RigidTerrain terrain(system);
    auto patch1_mat = chrono_types::make_shared<ChMaterialSurfaceSMC>();
    patch1_mat->SetFriction(friction);
    patch1_mat->SetRestitution(0.01f);
    auto patch1 = terrain.AddPatch(patch1_mat, ChCoordsys<>(ChVector<>(0, 0, 0), Q_from_AngX(slope)), 500, 500);
    std::cout << QUNIT << std::endl;
    patch1->SetColor(ChColor(1.0f, 1.0f, 1.0f));
    patch1->SetTexture(vehicle::GetDataFile("terrain/textures/tile4.jpg"), 200, 200);
    terrain.Initialize();

    ChFrame<> COMFrame = vehicle.GetCOMFrame();
    ChVector<> Pos = COMFrame.coord.pos;

    // Create the rigid body (this won't move, it is only for visualization tests)
    auto body = chrono_types::make_shared<ChBody>();
    body->SetBodyFixed(true);
    system->Add(body);
    auto orange_mat = chrono_types::make_shared<ChVisualMaterial>();
    orange_mat->SetDiffuseColor(ChColor(0.9f, 0.4f, 0.2f));
    auto cylinder = chrono_types::make_shared<ChCylinderShape>();
    cylinder->GetCylinderGeometry().p1 = ChVector<>(0.0,0.0, -0.5);
    cylinder->GetCylinderGeometry().p2 = ChVector<>(0.0,0.0, 0.5);
    cylinder->GetCylinderGeometry().rad = 3.0;
    cylinder->AddMaterial(orange_mat);
    body->AddVisualShape(cylinder, ChFrame<>(ChVector<>(vehicle.GetPointLocation(COMFrame.GetPos()).x(), vehicle.GetPointLocation(COMFrame.GetPos()).y(), vehicle.GetPointLocation(COMFrame.GetPos()).z()), QUNIT));


    auto body2 = chrono_types::make_shared<ChBody>();
    body2->SetBodyFixed(true);
    system->Add(body2);
    auto green_mat = chrono_types::make_shared<ChVisualMaterial>();
    green_mat->SetDiffuseColor(ChColor(0.0f, 1.0f, 0.0f));
    auto cylinder2 = chrono_types::make_shared<ChCylinderShape>();
    cylinder2->GetCylinderGeometry().p1 = ChVector<>(0.0, 0.0, -0.5);
    cylinder2->GetCylinderGeometry().p2 = ChVector<>(0.0, 0.0, 0.5);
    cylinder2->GetCylinderGeometry().rad = 3.6;
    cylinder2->AddMaterial(green_mat);
    body2->AddVisualShape(cylinder2, ChFrame<>(ChVector<>(goal[0], goal[1], 0), QUNIT));



    auto body3 = chrono_types::make_shared<ChBody>();
    body3->SetBodyFixed(true);
    system->Add(body3);
    auto black_mat = chrono_types::make_shared<ChVisualMaterial>();
    black_mat->SetDiffuseColor(ChColor(0.0f, 0.0f, 0.0f));
    //auto cylinder3 = chrono_types::make_shared<ChCylinderShape>();
    //cylinder3->GetCylinderGeometry().p1 = ChVector<>(0.0, 0.0, -9.0);
    //cylinder3->GetCylinderGeometry().p2 = ChVector<>(0.0, 0.0, 12.0);
    //cylinder3->GetCylinderGeometry().rad = 25.0;
    //cylinder3->AddMaterial(black_mat);
    //body3->AddVisualShape(cylinder3, ChFrame<>(ChVector<>(50.0, 50.0, -9.0 * sin(50.0 / 15.0) * sin(50.0 / 15.0)), QUNIT));

    auto cylinder4 = chrono_types::make_shared<ChCylinderShape>();
    cylinder4->GetCylinderGeometry().p1 = ChVector<>(0.0, 0.0, -9.0);
    cylinder4->GetCylinderGeometry().p2 = ChVector<>(0.0, 0.0, 12.0);
    cylinder4->GetCylinderGeometry().rad = 8.0;
    cylinder4->AddMaterial(black_mat);
    auto cylinder5 = chrono_types::make_shared<ChCylinderShape>();
    cylinder5->GetCylinderGeometry().p1 = ChVector<>(0.0, 0.0, 2.0);
    cylinder5->GetCylinderGeometry().p2 = ChVector<>(0.0, 0.0, -2.0);
    cylinder5->GetCylinderGeometry().rad = 2.5;
    cylinder5->AddMaterial(black_mat);


    body3->AddVisualShape(cylinder5, ChFrame<>(ChVector<>(50.0, -1.0, -1 * slope), QUNIT));
    body3->AddVisualShape(cylinder5, ChFrame<>(ChVector<>(70.0, 1.0, 1.0 * slope), QUNIT));
    body3->AddVisualShape(cylinder5, ChFrame<>(ChVector<>(90.0, -2.0, -1.0 * slope), QUNIT));
    //body3->AddVisualShape(cylinder5, ChFrame<>(ChVector<>(50.0, 90.0, -9.0 * sin(50.0 / 15.0) * sin(90.0 / 15.0)), QUNIT));
    //body3->AddVisualShape(cylinder5, ChFrame<>(ChVector<>(90.0, 50.0, -9.0 * sin(90.0 / 15.0) * sin(50.0 / 15.0)), QUNIT));


    std::string title = "Vehicle Control";
    double mycolor = 0.6;
    // Create the vehicle Irrlicht interface
    auto vis = chrono_types::make_shared<ChWheeledVehicleVisualSystemIrrlicht>();
    vis->SetWindowTitle(title);
    vis->SetWindowSize(1920, 1440);
    vis->SetChaseCamera(ChVector<>(-1.0, 0.0, 3.0), 6.0, 0.5);
    vis->Initialize();
    vis->AddLightDirectional();
    vis->AddSkyBox();
    vis->AddLogo();
    vis->AttachVehicle(&vehicle);


    // Create the interactive Irrlicht driver system
    //std::shared_ptr<ChDriver> driver;
    ChDriver driver = ChDriver(vehicle);
    driver.Initialize();

    DriverInputs driver_inputs = driver.GetInputs();
    int sim_count = 0;
    double mytime = sim_count * step_size;

    std::string outputfile_name = "MRZR_Controller_";
    std::cout << "output file name is: " << outputfile_name << std::endl;


    // initialize of file saving
    std::ofstream records;
    records.open("./MPC/" + outputfile_name + "records.csv", std::ios::out);

    std::ofstream tires;
    tires.open("./MPC/" + outputfile_name + "tires.csv", std::ios::out);

    std::ofstream wheel_pos;
    wheel_pos.open("./MPC/" + outputfile_name + "wheel_pos.csv", std::ios::out);


    vehicle.EnableRealtime(false);
    double start_time = 0;

    COMFrame = vehicle.GetCOMFrame();
    Pos = COMFrame.coord.pos;


    WheelState front_left_wheel = vehicle.GetWheel(0, VehicleSide::LEFT)->GetState();
    WheelState front_right_wheel = vehicle.GetWheel(0, VehicleSide::RIGHT)->GetState();
    WheelState rear_left_wheel = vehicle.GetWheel(1, VehicleSide::LEFT)->GetState();
    WheelState rear_right_wheel = vehicle.GetWheel(1, VehicleSide::RIGHT)->GetState();
    std::cout << vehicle.GetPointLocation(COMFrame.GetPos()).x() << std::endl;
    WheelState rear_wheel = vehicle.GetWheel(1, VehicleSide::LEFT)->GetState();
    std::cout << "COM Pos:  " << Pos << std::endl;
    WheelState front_wheel = vehicle.GetWheel(0, VehicleSide::LEFT)->GetState();


    std::cout << "Front Wheel Pos  " << front_wheel.pos << std::endl;

    std::cout << "Rear Wheel Pos  " << rear_wheel.pos << std::endl;
    vehicle.EnableRealtime(false);
    std::cout << vehicle.GetMass() << std::endl;
    std::cout << vehicle.GetInertia() << std::endl;

    bool control_flag = false;
    vehicle.EnableRealtime(false);

    //change here
    ChSpeedController speed_cont;
    speed_cont.SetGains(3.7, 0.2, 0.1);
    double targetSpeed = 4.0;
    double throttle_cmd = 0.0;

    if (img_output) {
        if (!filesystem::create_directory(filesystem::path(img_dir))) {
            std::cout << "Error creating directory " << img_dir << std::endl;
            return;
        }
    }
    int render_steps = (int)std::ceil(render_step_size / step_size);
    int step_number = 0;
    int render_frame = 0;

    double sa_list[3000] = { 0.0 };
    double ax_list[3000] = { 0.0 };
    double ux_list[3000] = { 4.0 };
    int control_idx = 0;
    double acc = 0;
    double ALPHA = 8.0;
    double GAMMA = 1.0;
    double acc_limit = 50;
    double acc_step = step_size * acc_limit;

    while (vis->Run()) {

        DriverInputs driver_inputs = driver.GetInputs();
        ChFrame<> COMFrame = vehicle.GetCOMFrame();
        double state[7];
        state[0] = vehicle.GetPointLocation(COMFrame.GetPos()).x();
        state[1] = vehicle.GetPointLocation(COMFrame.GetPos()).y();
        state[2] = vehicle.GetPointLocation(COMFrame.GetPos()).z();
        state[3] = vehicle.GetRot().e0();
        state[4] = vehicle.GetRot().e1();
        state[5] = vehicle.GetRot().e2();
        state[6] = vehicle.GetRot().e3();
        ChQuaternion<> q_temp(state[3], state[4], state[5], state[6]);
        ChVector<> velocity = vehicle.GetPointVelocity(vehicle.GetCOMFrame().GetPos());
        ChQuaternion<> q_trans = q_temp.GetInverse();
        velocity = q_trans.Rotate(velocity);
        sim_count = sim_count + 1;
        mytime = sim_count * step_size;

        ChVector<> velocityCOM = vehicle.GetPointVelocity(vehicle.GetCOMFrame().coord.pos);
        double Radius = vehicle.GetWheel(0, VehicleSide::LEFT)->GetRadius();
        velocity = q_trans.Rotate(velocityCOM);
        double ux = velocity.x();

        WheelState state0 = vehicle.GetWheel(0, VehicleSide::LEFT)->GetState();
        ChVector<> wheel_normal0 = state0.rot.GetYaxis();
        ChVector<> contact_point = ChVector<>(state0.pos.x(), state0.pos.y(), terrain.GetHeight(state0.pos));
        ChVector<> Z_dir0 = terrain.GetNormal(contact_point);
        ChVector<> X_dir0 = Vcross(wheel_normal0, Z_dir0);
        X_dir0.Normalize();
        ChVector<> Y_dir0 = Vcross(Z_dir0, X_dir0);
        ChMatrix33<> rot0;
        rot0.Set_A_axis(X_dir0, Y_dir0, Z_dir0);
        ChMatrix33<> rot0inv = rot0.transpose();
        double alpha_fl = vehicle.GetTire(0, VehicleSide::LEFT)->GetSlipAngle();
        double kappa_fl = vehicle.GetTire(0, VehicleSide::LEFT)->GetLongitudinalSlip();
        ChVector<> F_fl = rot0inv * (vehicle.GetTire(0, VehicleSide::LEFT)->ReportTireForce(&terrain).force);
        ChVector<> M_fl = rot0inv * (vehicle.GetTire(0, VehicleSide::LEFT)->ReportTireForce(&terrain).moment);

        WheelState state1 = vehicle.GetWheel(0, VehicleSide::RIGHT)->GetState();
        ChVector<> wheel_normal1 = state1.rot.GetYaxis();
        contact_point = ChVector<>(state1.pos.x(), state1.pos.y(), terrain.GetHeight(state1.pos));
        ChVector<> Z_dir1 = terrain.GetNormal(contact_point);
        ChVector<> X_dir1 = Vcross(wheel_normal1, Z_dir1);
        X_dir1.Normalize();
        ChVector<> Y_dir1 = Vcross(Z_dir1, X_dir1);
        ChMatrix33<> rot1;
        rot1.Set_A_axis(X_dir1, Y_dir1, Z_dir1);
        ChMatrix33<> rot1inv = rot1.transpose();
        double alpha_fr = vehicle.GetTire(0, VehicleSide::RIGHT)->GetSlipAngle();
        double kappa_fr = vehicle.GetTire(0, VehicleSide::RIGHT)->GetLongitudinalSlip();
        ChVector<> F_fr = rot1inv * (vehicle.GetTire(0, VehicleSide::RIGHT)->ReportTireForce(&terrain).force);
        ChVector<> M_fr = rot1inv * (vehicle.GetTire(0, VehicleSide::RIGHT)->ReportTireForce(&terrain).moment);

        WheelState state2 = vehicle.GetWheel(1, VehicleSide::LEFT)->GetState();
        ChVector<> wheel_normal2 = state2.rot.GetYaxis();
        contact_point = ChVector<>(state2.pos.x(), state2.pos.y(), terrain.GetHeight(state2.pos));
        ChVector<> Z_dir2 = terrain.GetNormal(contact_point);
        ChVector<> X_dir2 = Vcross(wheel_normal2, Z_dir2);
        X_dir2.Normalize();
        ChVector<> Y_dir2 = Vcross(Z_dir2, X_dir2);
        ChMatrix33<> rot2;
        rot2.Set_A_axis(X_dir2, Y_dir2, Z_dir2);
        ChMatrix33<> rot2inv = rot2.transpose();
        double alpha_rl = vehicle.GetTire(1, VehicleSide::LEFT)->GetSlipAngle();
        double kappa_rl = vehicle.GetTire(1, VehicleSide::LEFT)->GetLongitudinalSlip();
        ChVector<> F_rl = rot2inv * (vehicle.GetTire(1, VehicleSide::LEFT)->ReportTireForce(&terrain).force);
        ChVector<> M_rl = rot2inv * (vehicle.GetTire(1, VehicleSide::LEFT)->ReportTireForce(&terrain).moment);

        WheelState state3 = vehicle.GetWheel(1, VehicleSide::RIGHT)->GetState();
        ChVector<> wheel_normal3 = state3.rot.GetYaxis();
        contact_point = ChVector<>(state3.pos.x(), state3.pos.y(), terrain.GetHeight(state3.pos));
        ChVector<> Z_dir3 = terrain.GetNormal(contact_point);
        ChVector<> X_dir3 = Vcross(wheel_normal3, Z_dir3);
        X_dir3.Normalize();
        ChVector<> Y_dir3 = Vcross(Z_dir3, X_dir3);
        ChMatrix33<> rot3;
        rot3.Set_A_axis(X_dir3, Y_dir3, Z_dir3);
        ChMatrix33<> rot3inv = rot3.transpose();
        double alpha_rr = vehicle.GetTire(1, VehicleSide::RIGHT)->GetSlipAngle();
        double kappa_rr = vehicle.GetTire(1, VehicleSide::RIGHT)->GetLongitudinalSlip();
        ChVector<> F_rr = rot3inv * (vehicle.GetTire(1, VehicleSide::RIGHT)->ReportTireForce(&terrain).force);
        ChVector<> M_rr = rot3inv * (vehicle.GetTire(1, VehicleSide::RIGHT)->ReportTireForce(&terrain).moment);

        ChCoordsys<> tire_csys0(state0.pos, rot0.Get_A_quaternion());
        ChCoordsys<> tire_csys1(state1.pos, rot1.Get_A_quaternion());
        double steering = (Q_to_Euler123(tire_csys0.rot).z() - Q_to_Euler123(q_temp).z() + Q_to_Euler123(tire_csys1.rot).z() - Q_to_Euler123(q_temp).z()) / 2.;
        //if ((F_fl.z() == 0 || F_fr.z() == 0 || F_rl.z() == 0 || F_rr.z() == 0)) {
        //    std::cout << "Tire Lift Off:  " << mytime << std::endl;
        //}


        //std::cout << ux_list[control_idx]<<std::endl;
        double phi = Q_to_Euler123(q_temp).x();
        double theta = Q_to_Euler123(q_temp).y();
        double psi = Q_to_Euler123(q_temp).z();

        throttle_cmd = speed_cont.Advance(vehicle, ux_list[control_idx], step_size);
        double acc_ctl = 9.8 * (sin(phi)*sin(psi)-cos(phi)*cos(psi)*sin(theta)) + ALPHA * (ux_list[control_idx] - velocity.x()) + GAMMA * ax_list[control_idx];

        
        if (abs(acc - acc_ctl) > acc_step) {
            if (acc > acc_ctl) {
                acc = acc - acc_step;
            }
            else {
                acc = acc + acc_step;
            }
        }
        else {
            acc = acc_ctl;
        }
        throttle_cmd = std::max(std::min(acc / 3.2, 1.0), -1.0);

        //std::cout<< acc_ctl << "," << acc << ", " << throttle_cmd << std::endl;

        double braking_cmd = 0;
        if (throttle_cmd > 1.0) {
            throttle_cmd = 1.0;
            braking_cmd = 0;
        }
        else if (throttle_cmd < 0) {
            braking_cmd = -throttle_cmd;
            throttle_cmd = 0.0;
        }

        if (mytime <= WAIT_TIME) {
            driver_inputs.m_throttle = 0.0;
            driver_inputs.m_steering = 0.0;
            driver_inputs.m_braking = 1.0;
        }
        else {
            double str_ang = sa_list[control_idx];
            double steering_input = str_ang / 0.4554;
            driver_inputs.m_steering = steering_input;
            driver_inputs.m_throttle = throttle_cmd;
            driver_inputs.m_braking = braking_cmd;
        }
        control_idx = control_idx + 1;



        vis->BeginScene();
        vis->Render();
        vis->EndScene();

        if (img_output && step_number % render_steps == 0) {
            char filename[100];
            sprintf(filename, "%s/img_%03d.jpg", img_dir.c_str(), render_frame + 1);
            vis->WriteImageToFile(filename);
            render_frame++;
        }

        if ((ux <= 0 && mytime >= 5.0) || (  pow((state[0] -goal[0]),2) + pow((state[1] - goal[1]), 2) <= pow(3.6,2)   )) {
            vis->GetDevice()->closeDevice();
            std::ofstream terminationfile;
            terminationfile.open("E:/STUDY_UMICH/WIN2023/Local_Controller/stop.txt", std::ios::out);
            terminationfile.close();
            break;
        }

        step_number++;


        front_left_wheel = vehicle.GetWheel(0, VehicleSide::LEFT)->GetState();
        front_right_wheel = vehicle.GetWheel(0, VehicleSide::RIGHT)->GetState();
        rear_left_wheel = vehicle.GetWheel(1, VehicleSide::LEFT)->GetState();
        rear_right_wheel = vehicle.GetWheel(1, VehicleSide::RIGHT)->GetState();
        chrono::Vector compos = vehicle.GetPos() + vehicle.GetCOMFrame().GetPos();


        //// File IO Starts here
        records << mytime << "," << state[0] << "," << state[1] << "," << state[2] << "," << q_temp.e0() << "," << q_temp.e1() << "," << q_temp.e2() << "," << q_temp.e3() << "," << velocity.x() << "," << velocity.y() << "," << velocity.z() << "," << driver_inputs.m_steering * 0.4554 << "," << driver_inputs.m_throttle<<","<<ux_list[control_idx]<<","<<ax_list[control_idx] << std::endl;
        tires << alpha_fl << "," << kappa_fl << "," << F_fl.x() << "," << F_fl.y() << "," << F_fl.z() << "," << M_fl.x() << "," << M_fl.y() << "," << M_fl.z() << "," << alpha_fr << "," << kappa_fr << "," << F_fr.x() << "," << F_fr.y() << "," << F_fr.z() << "," << M_fr.x() << "," << M_fr.y() << "," << M_fr.z() << "," << alpha_rl << "," << kappa_rl << "," << F_rl.x() << "," << F_rl.y() << "," << F_rl.z() << "," << M_rl.x() << "," << M_rl.y() << "," << M_rl.z() << "," << alpha_rr << "," << kappa_rr << "," << F_rr.x() << "," << F_rr.y() << "," << F_rr.z() << "," << M_rr.x() << "," << M_rr.y() << "," << M_rr.z() << std::endl;

        if ((sim_count % 500) == 0) {
            std::ofstream cur_states;
            cur_states.open("E:/STUDY_UMICH/WIN2023/Local_Controller/" + outputfile_name + "cur_states.txt", std::ios::out);
            cur_states << mytime << "\t" << state[0] << "\t" << state[1]<<"\t" << Q_to_Euler123(q_temp).z() << "\t"  << velocity.y() <<"\t" << 0.0 <<"\t"<<velocity.x()<< "\t"<< driver_inputs.m_steering * 0.4554 << std::endl;
            cur_states.close();
            std::this_thread::sleep_for(std::chrono::milliseconds(1));
            control_idx = 0;
            while (true) {                
                std::ifstream send_signal;
                send_signal.open("cmd_sent.txt");
                if (send_signal) {
                    std::ifstream controls;
                    send_signal.close();
                    controls.open("MPC_cmd.txt");
                    std::cout << "The current time is:"<< mytime <<"; " << "MPC file exists" << std::endl;
                    for (int readidx = 0; readidx < 500; readidx++) {
                        controls >> sa_list[readidx];
                        controls >> ax_list[readidx];
                        controls >> ux_list[readidx];
                        //std::cout << sa_list[readidx] << " " << ax_list[readidx] << " " << ux_list[readidx] << std::endl;
                    }
                    controls.close();
                    std::this_thread::sleep_for(std::chrono::milliseconds(1));
                    std::remove("MPC_cmd.txt");
                    std::remove("cmd_sent.txt");
                    break;
                }
                else {
                    std::this_thread::sleep_for(std::chrono::milliseconds(1));
                }
            }
        }





        // Update modules (process inputs from other modules)
        double time = vehicle.GetSystem()->GetChTime();
        driver.Synchronize(time);
        vehicle.Synchronize(time, driver_inputs, terrain);
        terrain.Synchronize(time);
        vis->Synchronize(time, driver_inputs);

        // Advance simulation for one timestep for all modules
        driver.Advance(step_size);
        vehicle.Advance(step_size);
        terrain.Advance(step_size);
        vis->Advance(step_size);
    }
    records.close();
    tires.close();
    wheel_pos.close();
    vis->GetDevice()->closeDevice();
}



int main(int argc, char* argv[]) {
    run_sim();
    return 0;
}