// Builds a graph with errors, displays it and saves it
// as image. First, include some header files (within,
// CINT these will be ignored).

#include "TArrow.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TROOT.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TH1I.h"
#include "TObject.h"
#include <assert.h>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iostream>
#include <map>
#include <sstream>
#include <vector>

// Returns the column names line, i.e. the last line of the header.
// After calling this method the stream is ready to read data.
std::string read_header(std::istream& file_in, int header_lines = 1);
// Returns the column names line, i.e. the last line of the header.
std::string read_header(std::string const& file_name, int header_lines = 1);
// Returns the column names line, i.e. the last line of the header.
std::string read_header(std::filesystem::path const& file_path,
                        int header_lines = 1);

// Returns a vector containing the elements separated by 'sep' in 'line'
std::vector<std::string> split(std::string const& line, char const sep = ' ');

// Returns the value of the 'index' element in a vector of strings representing
// numbers
// This is where conversion from string to value happens.
template<class T>
T get_line_element(std::vector<std::string> const& line_vec, int index);

// Returns a vector containg values read from the 'index' column of
// 'file_in_name'
template<class T>
std::vector<T> get_column(std::string const& file_in_name, int index);
template<class T>
std::vector<T> get_column(std::string const& file_in_name,
                          std::string const& name);

struct PlotAttributes
{
  std::string plot_function{""};
  std::string draw_options{""};
};

void progress(std::ostream& out = std::cout)
{
  std::vector<char> progress_spin = {'|', '/', '-', '\\'};
  for (auto& val : progress_spin) {
    out << '\r' << val << std::flush;
  }
}
// TArgs must be types for which operator<< is overloaded
template<class... TArgs>
void flush_message(std::ostream& out, TArgs const&... out_args)
{
  (out << ... << out_args) << std::flush;
}

// Plot using a map containing vectors of data.
template<class T>
void plotxt(std::map<std::string, std::vector<T>> axes,
            std::vector<std::string> labels, PlotAttributes const& plot_attr)
{
  assert(!axes.empty() && !labels.empty() && "No data to plot.");

  flush_message(std::cout, "\nPlotting...");

  if (plot_attr.plot_function == "TGraph") {
    TGraph* plot = new TGraph(axes[labels[0]].size(), axes[labels[0]].data(),
                              axes[labels[1]].data());

    plot->Draw(plot_attr.draw_options.c_str());

  } else if (plot_attr.plot_function == "TH1I") {
    auto const [min_it, max_it] =
        std::minmax_element(axes[labels[0]].begin(), axes[labels[0]].end());

    std::string title{static_cast<std::string>("Title;") + labels[0]
                      + ";Counts"};
    TH1I* plot = new TH1I("hist", title.c_str(), axes[labels[0]].size(),
                          *min_it, *max_it);
    for (int const val : axes[labels[0]]) {
      plot->Fill(val);
    }

    plot->Draw(plot_attr.draw_options.c_str());

  }
  /* //Add plot functions
  else if (plot_attr.plot_function == ??) {
    // build plot object
    plot->Draw(plot_attr.draw_options);
  }*/
  else {
    throw std::invalid_argument("\nPlot " + plot_attr.plot_function
                                + " is not supported");
  }
}

// Plot passing column indexes as parameters
template<class T>
void plotxt(std::filesystem::path data_file,
            std::vector<int> const& column_indexes,
            PlotAttributes const& plot_attr)
{
  std::map<std::string, std::vector<T>> axes;
  std::vector<std::string> labels;

  auto column_names = split(read_header(data_file));
  for (int const index : column_indexes) {
    assert(index < column_names.size()
           && "Index is greater than available columns");
    std::string label = column_names[index];
    labels.push_back(label);
    axes[label] = get_column<T>(data_file, index);
  }

  return plotxt<T>(axes, labels, plot_attr);
}

// Required to pass brace-list initializer as argument
template<class T>
void plotxt(std::filesystem::path data_file,
            std::initializer_list<int> const& column_indexes,
            PlotAttributes const& plot_attr)
{
  std::vector<int> indexes_vec = column_indexes;
  return plotxt<T>(data_file, indexes_vec, plot_attr);
}

// Plot using column names as parameters
template<class T>
void plotxt(std::filesystem::path data_file,
            std::vector<std::string> const& column_names,
            PlotAttributes const& plot_attr)
{
  std::map<std::string, std::vector<T>> axes;
  std::vector<std::string> labels = column_names;
  for (std::string const& label : column_names) {
    axes[label] = get_column<T>(data_file, label);
  }

  return plotxt<T>(axes, labels, plot_attr);
}

// Required to pass brace-list initializer as argument
template<class T>
void plotxt(std::filesystem::path data_file,
            std::initializer_list<std::string> const& column_indexes,
            PlotAttributes const& plot_attr)
{
  std::vector<std::string> names_vec = column_indexes;
  return plotxt<T>(data_file, names_vec, plot_attr);
}

#ifndef __CINT__
#include "TApplication.h"
#include <boost/program_options.hpp>
#include <algorithm>
#include <iomanip>
#include <thread>
namespace po = boost::program_options;

// TArgs must be types for which operator<< is overloaded
template<class... TArgs>
void exit_message(std::ostream& out, TArgs const&... out_args)
{
  (out << ... << out_args) << std::endl;
  std::exit(-1);
}

void StandaloneApplication(int argc, char** argv)
{
  std::filesystem::path data_file;
  bool list_columns{false};
  std::vector<std::string> data_columns;
  PlotAttributes plot_attr;

  std::string usage =
  // clang-format off
    static_cast<std::string>("Usage:") // Cast needed to use the operator+
    + "\n  ./pltxt -h" + "\n  ./pltxt -d <data_file> -l"
    + "\n  ./pltxt -d <data_file> -c COLUMN_INDEX | COLUMN_NAME [COLUMN_INDEX | COLUMN_NAME...] -p PLOT_FUNCTION [-o DRAW_OPTIONS]";
  // clang-format on

  po::options_description desc{"Options"};
  desc.add_options()
    // clang-format off
  ("help,h", "Print help message.")
  ("data-file,d", po::value<std::filesystem::path>(&data_file)->value_name("<data_file>")->required(), "The file containing data arranged in columns.")
  ("list-column-names,l", po::bool_switch(&list_columns)->default_value(false), "Show column names. Read them from the first line before data in <data_file>.")
  ("data-columns,c", po::value<std::vector<std::string>>(&data_columns)->multitoken()->value_name("COLUMN_NAME | COLUMN_INDEX"), "The indexes or names of columns containing data to plot.")
  ("plot-function,p", po::value<std::string>(&plot_attr.plot_function)->value_name("PLOT_FUNCTION"), "The ROOT function for plotting data. See ROOT documentation to know more.")
  ("draw-options,o", po::value<std::string>(&plot_attr.draw_options)->value_name("DRAW_OPTIONS")->default_value(""), "Options to pass to ROOT Draw function. See ROOT documentation to know more.");
  // clang-format on

  try {
    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).options(desc).run(), vm);
    if (vm.count("help")) {
      exit_message(std::cout, usage, "\n\n", desc);
    }

    po::notify(vm);
  } catch (po::error const& e) {
    exit_message(std::cerr, e.what(), '\n', usage, '\n');
  } catch (std::exception const& e) {
    exit_message(std::cerr, e.what());
  } catch (...) {
    exit_message(std::cerr, "Unknown exception. Quitting...");
  }


  if (list_columns) {
    exit_message(std::cout, read_header(data_file));
  }

  if (data_columns.empty()) {
    exit_message(std::cout, usage);
  }

  flush_message(std::clog, "\nReading data...");
  try {
    std::vector<int> index_vec;
    for (std::string const& column : data_columns) {
      index_vec.push_back(std::stoi(column));
    }

    myplot<int>(data_file, index_vec, plot_attr);
  } catch (std::invalid_argument& ia) {
    myplot<int>(data_file, data_columns, plot_attr);
  }

    //macro1();
}

int main(int argc, char** argv)
{
  // Do not let ROOT receive command line arguments.
  TApplication app("ROOT Application", 0, {});
  StandaloneApplication(argc, argv);
  app.Run();
  return 0;
}
#endif


std::string read_header(std::istream& file_in, int header_lines)
{
  std::string line_str;
  // Label line is the last line of the header
  do
    std::getline(file_in, line_str);
  while (header_lines > 1);
  return line_str;
}
std::string read_header(std::string const& file_name, int header_lines)
{
  std::ifstream file_in(file_name);
  return read_header(file_in, header_lines);
}
std::string read_header(std::filesystem::path const& file_name,
                        int header_lines)
{
  std::ifstream file_in(file_name);
  return read_header(file_in, header_lines);
}

std::vector<std::string> split(std::string const& line, char const sep)
{
  std::istringstream line_in(line);
  std::vector<std::string> line_vector{};
  std::string line_element{};

  while (std::getline(line_in, line_element, sep))
    line_vector.push_back(line_element);
  return line_vector;
}

// This is where conversion from string to value happens.
template<class T>
T get_line_element(std::vector<std::string> const& line_vec, int index);

template<>
int get_line_element<int>(std::vector<std::string> const& line_vec, int index)
{
  assert(index < line_vec.size() && "Index is out of range.");
  return std::stoi(line_vec.at(index));
}

template<>
double get_line_element<double>(std::vector<std::string> const& line_vec,
                                int index)
{
  assert(index < line_vec.size() && "Index is out of range.");
  return std::stod(line_vec.at(index));
}

template<class T>
std::vector<T> get_column(std::string const& file_in_name, int index)
{
  std::ifstream file_in(file_in_name);
  read_header(file_in);

  std::vector<T> column{};
  std::string line_str;
  while (std::getline(file_in, line_str)) {
    T val = get_line_element<T>(split(line_str), index);
    column.push_back(val);
  }
  return column;
}

template<class T>
std::vector<T> get_column(std::string const& file_in_name,
                          std::string const& name)
{
  auto labels   = split(read_header(file_in_name));
  auto label_it = std::find(labels.begin(), labels.end(), name);
  if (label_it == labels.end())
    throw std::invalid_argument("Column name not found in file" + file_in_name);

  int index = std::distance(labels.begin(), label_it);

  return get_column<T>(file_in_name, index);
}
