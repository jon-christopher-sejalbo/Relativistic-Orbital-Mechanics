#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>

// Defining Constants [SI UNITS]
const double c = 2.99e8; // Speed of light
const double M_BH = 3*1.988e31; // Black hole mass; currently 30 solar masses
const double G = 6.67e-11; // Gravitational constant
const double M_BH_GEOM = (G*M_BH)/(c*c); // Geometric mass of the black hole
const double e = 1.05; // Relativistic energy

// Quadratic Components for r_c
const double A = (1-e*e);
const double B = (3*M_BH_GEOM*e*e-4*M_BH_GEOM);
const double C = 4*M_BH_GEOM*M_BH_GEOM;

const double r_c = (-B-pow(B*B-4*A*C,0.5))/(2*A); // Critical radius
const double l = pow((r_c*M_BH_GEOM)/(1-(3*M_BH_GEOM)/r_c),0.5)+(1e-5*M_BH_GEOM); // Angular momentum
const double b = (l/pow((e*e-1),0.5)); // Impact parameter
const double beta = pow((e*e-1),0.5); // Fraction of the speed of light
const double lorentz_factor = pow(1-beta*beta,-0.5); // Lorentz factor
const double r_init = 50*r_c; // Initial radius

// Time stepping constants
const double t_tot = 1e10;
const double tol = 1e-12;
const double scaling_factor = 1.0001; // Scaling factor for adaptive integrator
const double max_dtau = 0.0001; // Max value for dtau
const double min_dtau = 1e-30; // Minimum value for dtau
const double shrink_factor = 0.5; // Shrink factor for adaptive integrator

// Stepping variables
const int step = 1;

struct State {
    double x_1, y_1, x_2, y_2;
};  

struct Vector4 {
    double coordinates[4];

    // Default constructor
    Vector4() {
        for(int i=0;i<4;i++) {
            this->coordinates[i] = 0;
        }
    }

    // Parameter constructor
    Vector4(double t, double x, double y, double z) {
        this->coordinates[0] = t;
        this->coordinates[1] = x;
        this->coordinates[2] = y;
        this->coordinates[3] = z;
    }

    // Copy constructor
    Vector4(const Vector4 &other) {
        for(int i=0;i<4;i++) {
            this->coordinates[i] = other.coordinates[i];
        }
    }

    // Prints out the vector in <...> notation
    void print_Vector(std::string title) {
        std::cout<<title<<": <";
        for(int i=0;i<4;i++) {
            std::cout<<this->coordinates[i];
            if(i!=3) {
                std::cout<<",";
            }
        }
        std::cout<<">"<<std::endl;
    }
};


struct Orbital_Mechanics {
    // Updates angular velocity
    double update_angular_velocity(double r, double l) {
        return l/(r*r);
    }

    // Updates radial velocity
    double update_radial_velocity(double r, double rel_energy_mass, double l, bool is_ingoing) {
        double inner_term = (rel_energy_mass)*(rel_energy_mass) - 1 + (2*M_BH_GEOM*r*r-l*l*r+2*M_BH_GEOM*l*l)/(r*r*r);
        if(inner_term<0) {
            std::cout<<r<<std::endl;
            inner_term=0;
        }
        double radial_vel = pow((inner_term),0.5);
        return is_ingoing ? -radial_vel : radial_vel;
    }
};

// Polar-to-Cartesian coordinate converter
Vector4 polar_Cartesian(const Vector4 &polar) {
    Vector4 resultant = Vector4();

    resultant.coordinates[0] = polar.coordinates[0];
    resultant.coordinates[1] = polar.coordinates[1]*cos(polar.coordinates[2]);
    resultant.coordinates[2] = polar.coordinates[1]*sin(polar.coordinates[2]);
    resultant.coordinates[3] = 0; // Equitorial plane condition 

    return resultant;
}

int main() {
    bool is_ingoing = true;

    // Data vector
    std::vector<State> states;
    states.reserve(10000000);

    // Stepping variables
    double dtau = 0.0000001;
    int k = 0;

    // Coordinate time variable
    double t_coord = 0;
    double t_proper = 0;

    // Initializing orbital mechanics
    Orbital_Mechanics mechanics = Orbital_Mechanics();

    // Rocket vectors
    Vector4 rocket_position(0,r_init,0,0);
    Vector4 rocket_velocity(0,mechanics.update_radial_velocity(r_init, e, l, is_ingoing),mechanics.update_angular_velocity(r_init, l),0);

    bool found_anom = false;

    while(t_proper<t_tot) {
        // Conversion to coordinate time
        double dt = ((1/(((1-((2*M_BH_GEOM)/(rocket_position.coordinates[1]))))/(lorentz_factor)))*dtau);

        // Data entries
        if(k%step==0) {
            Vector4 cartesian = polar_Cartesian(rocket_position);
            states.push_back({cartesian.coordinates[1],cartesian.coordinates[2],0,0});
        }

        // Event horizon capture check
        if(rocket_position.coordinates[1]<=(2*M_BH_GEOM)) {
            std::cout<<"Captured by event horizon, goodbye spaceship!"<<std::endl;
            break;
        }

        // Adaptive time-stepping implementation
        double init_pos = rocket_position.coordinates[1];
        double init_vel = rocket_velocity.coordinates[1];
        double pred_pos = init_pos + init_vel*dtau;
        double curr_vel = init_vel;
        double curr_pos = init_pos;
        double init_phi = rocket_position.coordinates[2];
        double pred_phi = init_phi + mechanics.update_angular_velocity(init_pos, l)*dtau;
        double curr_phi = init_phi;

        for(int i=1;i<3;i++) {
            curr_vel  = mechanics.update_radial_velocity(curr_pos, e, l, is_ingoing);
            curr_pos += curr_vel*0.5*dtau;
            curr_phi += mechanics.update_angular_velocity(curr_pos, l)*0.5*dtau;
        }

        if(abs(curr_pos-pred_pos)<tol&&abs(curr_phi-pred_phi)<tol) {
            rocket_velocity.coordinates[1] = curr_vel;
            rocket_position.coordinates[1] = curr_pos;
            if(curr_pos>=init_pos&&is_ingoing) {
                std::cout<<"FLAG!"<<std::endl;
                std::cout<<"FLAGGED RADIUS: "<<curr_pos<<std::endl;
                std::cout<<"FLAGGED VELOCITY: "<<curr_vel<<std::endl;
                is_ingoing = false;
                rocket_position.coordinates[1]+=1e-5*M_BH_GEOM;
            }
            rocket_position.coordinates[2] = curr_phi;
            t_coord +=dt;
            t_proper+=dtau;
            dtau*=scaling_factor;
        } else {
            if(dtau<min_dtau) {
                dtau=min_dtau;
                rocket_position.coordinates[1] = curr_pos;
                rocket_position.coordinates[2] = curr_phi;
                t_proper+=dtau;
            }
            std::cout<<"Rejected Step!"<<std::endl;
            dtau*=shrink_factor;
        }

        k++;
    }

    // Proper-time coordinate-time comparison
    std::cout<<"Coordinate time to an observer at r->INF: "<<t_coord<<"s"<<std::endl;
    std::cout<<"Your wristwatch/proper time: "<<t_proper<<"s"<<std::endl;
    std::cout<<"Time dialation: "<<((t_coord-t_proper)/(t_proper))*100<<"%"<<std::endl;
    std::cout<<"Cool, right?"<<std::endl;

    std::ofstream file("deep_space_whirling.csv");
    for(const auto &pt : states) {
        file<<pt.x_1<<","<<pt.y_1<<","<<pt.x_2<<","<<pt.y_2<<'\n';
    }

    return 0;
}