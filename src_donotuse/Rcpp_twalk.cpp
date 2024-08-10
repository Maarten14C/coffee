#include <Rcpp.h>
using namespace Rcpp;
using namespace std;


double IntProd(NumericVector x) {
  return sum(x * x);
}

double DotProd(NumericVector x, NumericVector y) {
  return sum(x * y);
}


// h1 function for the "traverse" kernel 
List Simh1(int dim, double pphi, NumericVector x, NumericVector xp, double beta) {
  GetRNGstate();
  NumericVector rt(dim);
  NumericVector ran_vec = Rcpp::runif(dim);
  LogicalVector phi = ran_vec < pphi;
  PutRNGstate(); // Ensure RNG state is save
  for(int i = 0; i < dim; i++) {
    if(phi[i])
      rt[i] = xp[i] + beta * (xp[i] - x[i]);
    else
      rt[i] = x[i];
  }
  double nphi = sum(phi);
  return List::create(Named("rt") = rt, Named("nphi") = nphi );
}

//Simulation of the beta parameter for kernel h1          
double Simfbeta(double at) {
  GetRNGstate();
  double rand_number = R::runif(0, 1);
  PutRNGstate(); // Ensure RNG state is saved
  GetRNGstate();
  double rand_number2 = R::runif(0, 1);
  PutRNGstate(); // Ensure RNG state is saved
  if(rand_number < (at - 1) / (2 * at))
    return exp(log(rand_number2) / (at + 1));
  else
    return exp(log(rand_number2) / (1 - at));
}

//h function for the "walk" kernel 
List Simh2(int dim, double pphi, double aw, NumericVector x, NumericVector xp) {
  GetRNGstate();
  NumericVector u = Rcpp::runif(dim);
  PutRNGstate(); // Ensure RNG state is saved
  GetRNGstate();
  NumericVector ran_vec = Rcpp::runif(dim);
  PutRNGstate(); // Ensure RNG state is saved
  LogicalVector phi = ran_vec < pphi;
  // Convert LogicalVector to NumericVector
  NumericVector phi_numeric = as<NumericVector>(phi);
  NumericVector z = ((aw / (1 + aw)) * (aw * pow(u, 2) + 2 * u - 1)) * phi_numeric;
  // z = z * phi_numeric;
  NumericVector rt = x + (x - xp) * z;
  double nphi = sum(phi);
  return List::create(Named("rt") = rt, Named("nphi") = nphi);
}

//h function for the "blow" kernel 
List Simh3(int dim, double pphi, NumericVector x, NumericVector xp) {
  GetRNGstate();
  NumericVector ran_vec = Rcpp::runif(dim);
  LogicalVector phi = ran_vec < pphi;
  PutRNGstate(); // Ensure RNG state is saved
  GetRNGstate();
  NumericVector ran_norm = Rcpp::rnorm(dim);
  PutRNGstate(); // Ensure RNG state is saved
  // Convert LogicalVector to NumericVector
  NumericVector phi_numeric = as<NumericVector>(phi);
  NumericVector diff = abs(xp - x);
  double sigma = max(diff * phi_numeric);
  NumericVector rt = xp * phi_numeric + sigma * ran_norm * phi_numeric + x * (1 - phi_numeric);
  double nphi = sum(phi);
  return List::create(Named("rt") = rt, Named("nphi") = nphi, Named("phi") = phi);
}

//-log of g3, the density of h_b
double G3U(double nphi, LogicalVector phi, NumericVector h, NumericVector x, NumericVector xp) {
  NumericVector diff = abs(xp - x);
  // Convert LogicalVector to NumericVector
  NumericVector phi_numeric = as<NumericVector>(phi);
  double sigma = max(diff * phi_numeric);
  if(nphi > 0) {
    return (nphi / 2.0) * log(2.0 * M_PI) + nphi * log(sigma) + 0.5 * IntProd(h - xp) / (sigma * sigma);
  } else {
    return 0.0;
  }
}


//h function for the "hop" kernel 
List Simh4(int dim, double pphi, NumericVector x, NumericVector xp) {
  GetRNGstate();
  NumericVector ran_vec = Rcpp::runif(dim);
  LogicalVector phi = ran_vec < pphi;
  PutRNGstate(); // Ensure RNG state is saved
  GetRNGstate();
  double ran_norm = R::rnorm(0, 1);
  PutRNGstate(); // Ensure RNG state is saved
  // Convert LogicalVector to NumericVector
  NumericVector phi_numeric = as<NumericVector>(phi);
  NumericVector diff = abs(xp - x);
  double sigma = max(diff * phi_numeric) / 3.0;
  NumericVector rt(dim);
  for(int i = 0; i < dim; i++) {
    if(phi[i]) {
      rt[i] = x[i] + sigma * ran_norm;
    } else {
      rt[i] = x[i];
    }
  }
  int nphi = sum(phi);
  return List::create(Named("rt") = rt, Named("nphi") = nphi, Named("phi") = phi);
}

const double log2pi = log(2 * M_PI);
const double log3 = log(3);

//-log of g4, the density of h_h.
double G4U(double nphi, LogicalVector phi, NumericVector h, NumericVector x, NumericVector xp) {
  NumericVector diff = abs(xp - x);
  // Convert LogicalVector to NumericVector
  NumericVector phi_numeric = as<NumericVector>(phi);
  double sigma = max(diff * phi_numeric) / 3.0;
  if(nphi > 0) {
    return (nphi / 2.0) * log2pi - nphi * log3 + nphi * log(sigma) + 0.5 * 9 * IntProd(h - x) / (sigma * sigma);
  } else {
    return 0.0;
  }
}

//////////////////////////
// One move

// Assume the functions Simh1, Simh2, Simh3, Simh4, G3U, G4U, and Simfbeta are already defined

// [[Rcpp::export]]
List OneMove(int dim, Function Obj, Function Supp, NumericVector x, double U, NumericVector xp, double Up, 
             double at=6, double aw=1.5, double pphi=0.5, double F1=0.4918, double F2=0.9836, double F3=0.9918) {
  GetRNGstate();
  double ker = R::runif(0, 1);
  PutRNGstate(); // Ensure RNG state is saved
  
  NumericVector y(dim);
  NumericVector yp(dim);
  double propU;
  double propUp;
  double A;
  double nphi;
  LogicalVector phi;
  GetRNGstate();
  double r_value = R::runif(0, 1) ;
  PutRNGstate(); // Ensure RNG state is saved
  
  if(ker < F1) { // Kernel h1: traverse
    double beta = Simfbeta(at);
    
    if(r_value < 0.5) {
      List tmp = Simh1(dim, pphi, xp, x, beta);
      yp = tmp["rt"];
      nphi = tmp["nphi"];
      y = clone(x);
      propU = U;
      if(as<bool>(Supp(yp))) {
        propUp = as<double>(Obj(yp));
        if(nphi == 0) {
          A = 1;
        }
        else{
          A = exp((U - propU) + (Up - propUp) + (nphi - 2) * log(beta));
        }
      } else {
        propUp = NA_REAL;
        A = 0;
      }
    } else {
      List tmp = Simh1(dim, pphi, x, xp, beta);
      y = tmp["rt"];
      nphi = tmp["nphi"];
      yp = clone(xp);
      propUp = Up;
      if(as<bool>(Supp(y))) {
        propU = as<double>(Obj(y));
        if(nphi == 0) {
          A = 1;
          }
        else {
          A = exp((U - propU) + (Up - propUp) + (nphi - 2) * log(beta));
          }
      } else {
        propU = NA_REAL;
        A = 0;
      }
    }
  } else if(ker < F2) { // Kernel h2: walk

    if(r_value < 0.5) {
      List tmp = Simh2(dim, pphi, aw, xp, x);
      yp = tmp["rt"];
      nphi = tmp["nphi"];
      y = clone(x);
      propU = U;
      if(as<bool>(Supp(yp)) && is_true(all(abs(yp - y) > 0))) {
        propUp = as<double>(Obj(yp));
        A = exp((U - propU) + (Up - propUp));
      } else {
        propUp = NA_REAL;
        A = 0;
      }
    } else {
      List tmp = Simh2(dim, pphi, aw, x, xp);
      y = tmp["rt"];
      nphi = tmp["nphi"];
      yp = clone(xp);
      propUp = Up;
      if(as<bool>(Supp(y)) && is_true(all(abs(yp - y) > 0))) {
        propU = as<double>(Obj(y));
        A = exp((U - propU) + (Up - propUp));
      } else {
        propU = NA_REAL;
        A = 0;
      }
    }
  } else if(ker < F3) { // Kernel h3: blow
    if(r_value < 0.5) {
      List tmp = Simh3(dim, pphi, xp, x);
      yp = tmp["rt"];
      nphi = tmp["nphi"];
      phi = tmp["phi"];
      y = clone(x);
      propU = U;
      if(as<bool>(Supp(yp)) && is_true(all(yp != x))) {
        propUp = as<double>(Obj(yp));
        double W1 = G3U(nphi, phi, yp, xp, x);
        double W2 = G3U(nphi, phi, xp, yp, x);
        A = exp((U - propU) + (Up - propUp) + (W1 - W2));
      } else {
        propUp = NA_REAL;
        A = 0;
      }
    } else {
      List tmp = Simh3(dim, pphi, x, xp);
      y = tmp["rt"];
      nphi = tmp["nphi"];
      phi = tmp["phi"];
      yp = clone(xp);
      propUp = Up;
      if(as<bool>(Supp(y)) && is_true(all(y != xp))) {
        propU = as<double>(Obj(y));
        double W1 = G3U(nphi, phi, y, x, xp);
        double W2 = G3U(nphi, phi, x, y, xp);
        A = exp((U - propU) + (Up - propUp) + (W1 - W2));
      } else {
        propU = NA_REAL;
        A = 0;
      }
    }
  } else { // Kernel h4: hop
    if(r_value < 0.5) {
      List tmp = Simh4(dim, pphi, xp, x);
      yp = tmp["rt"];
      nphi = tmp["nphi"];
      phi = tmp["phi"];
      y = clone(x);
      propU = U;
      if(as<bool>(Supp(yp)) && is_true(all(yp != x))) {
        propUp = as<double>(Obj(yp));
        double W1 = G4U(nphi, phi, yp, xp, x);
        double W2 = G4U(nphi, phi, xp, yp, x);
        A = exp((U - propU) + (Up - propUp) + (W1 - W2));
      } else {
        propUp = NA_REAL;
        A = 0;
      }
    } else {
      List tmp = Simh4(dim, pphi, x, xp);
      y = tmp["rt"];
      nphi = tmp["nphi"];
      phi = tmp["phi"];
      yp = clone(xp);
      propUp = Up;
      if(as<bool>(Supp(y)) && is_true(all(y != xp))) {
        propU = as<double>(Obj(y));
        double W1 = G4U(nphi, phi, y, x, xp);
        double W2 = G4U(nphi, phi, x, y, xp);
        A = exp((U - propU) + (Up - propUp) + (W1 - W2));
      } else {
        propU = NA_REAL;
        A = 0;
      }
    }
  }
  // Rcout << "val aletorio: " << r_value << std::endl;
  // Rcout << "supp supp fun: " << as<bool>(Supp(y)) << std::endl;
  // Rcout << "y: " << y << std::endl;


  return List::create(Named("y") = y, Named("yp") = yp, Named("A") = A, Named("propU") = propU, Named("propUp") = propUp);
}
///////////////////////
// Runtwalk function

// Assume the functions Simh1, Simh2, Simh3, Simh4, G3U, G4U, Simfbeta, and OneMove are already defined

// [[Rcpp::export]]
List Runtwalk(int dim, Function Obj, Function Supp, NumericVector x0, NumericVector xp0, int niter, int burnin, int thin,
              double at=6, double aw=1.5, double pphi=0.5, double F1=0.4918, double F2=0.9836, double F3=0.9918) {
  NumericVector x = clone(x0);
  NumericVector xp = clone(xp0);
  double U = as<double>(Obj(x));
  double Up = as<double>(Obj(xp));
  
  NumericMatrix samples(niter, dim);
  NumericMatrix samplesp(niter, dim);
  NumericVector Us(niter);
  NumericVector Ups(niter);
  
  // Initialize the progress bar
  int total_iter = burnin +  niter * thin;
  int progress = 0;
  
  // Print the message
  Rcout << "Total iterations: " << total_iter << "\n";
  Rcout << "Burn-in period: " << burnin << "\n";
  Rcout << "Saving every " << thin << " iterations\n";
  Rcout << "Number of iterations to be saved: " << niter << "\n";
  
  // Function to print the progress bar
  auto print_progress = [](int progress, int total) {
    int width = 50; // Width of the progress bar
    int pos = width * progress / total;
    Rcout << "[";
    for (int i = 0; i < width; ++i) {
      if (i < pos) Rcout << "=";
      else if (i == pos) Rcout << ">";
      else Rcout << " ";
    }
    Rcout << "] " << int(progress / (double)total * 100.0) << " %\r";
    Rcout.flush();
  };
  
  for(int burn = 0; burn < burnin; burn++) {
    List move = OneMove(dim, Obj, Supp, x, U, xp, Up, at, aw, pphi, F1, F2, F3);
    Vector move_y = as<NumericVector>(move["y"]);
    Vector move_yp = as<NumericVector>(move["yp"]);
    double move_U = as<double>(move["propU"]);
    double move_Up = as<double>(move["propUp"]);
    double move_A = as<double>(move["A"]);
    
    GetRNGstate();
    double r_value = R::runif(0, 1) ;
    PutRNGstate(); // Ensure RNG state is saved
    if(r_value<  move_A ) {
      x = move_y;
      xp = move_yp;
      U = move_U;
      Up = move_Up;
    }
    // Update and print the progress bar
    progress = progress + 1;
    print_progress(progress, total_iter);
  }
  
  // Announce completion of burn-in phase
  Rcout << std::endl << "Burn-in period completed." << std::endl;
  
  for(int iter = 0; iter < niter; ) {
    for(int thinIter = 0; thinIter < thin; thinIter++) {
      List move = OneMove(dim, Obj, Supp, x, U, xp, Up, at, aw, pphi, F1, F2, F3);
      Vector move_y = as<NumericVector>(move["y"]);
      Vector move_yp = as<NumericVector>(move["yp"]);
      double move_U = as<double>(move["propU"]);
      double move_Up = as<double>(move["propUp"]);
      double move_A = as<double>(move["A"]);

      GetRNGstate();
      double r_value = R::runif(0, 1) ;
      PutRNGstate(); // Ensure RNG state is saved
      if(r_value < move_A ) {
        x = move_y;
        xp = move_yp;
        U = move_U;
        Up = move_Up;
      }
      // Update and print the progress bar
      progress = progress + 1;
      print_progress(progress, total_iter);

    }
    
    for(int j = 0; j < dim; j++) {
      samples(iter, j) = x[j];
      samplesp(iter, j) = xp[j];
    }
    Us(iter) = U;
    Ups(iter) = Up;

    iter++;

    
    // Debugging
    // Print values at specific intervals
    // if (iter % (niter / 10) == 0) {
    //   Rcout << "Iteration: " << iter << std::endl;
    //   Rcout << "x: " << x << std::endl;
    //   Rcout << "xp: " << xp << std::endl;
    //   Rcout << "U: " << U << std::endl;
    //   Rcout << "Up: " << Up << std::endl;
    // }

  }
  
  // Finalize the progress bar
  Rcout << std::endl;
  
  
  return List::create(Named("samples") = samples, Named("samplesp") = samplesp, Named("Us") = Us, Named("Ups") = Ups);
}











