#include "graphio.hpp"
#include <boost/algorithm/string.hpp>
#include <iostream>
#include <map>
#include <tinyxml2.h>

namespace pds {

struct ParseError : std::exception {
private:
  std::string reason;

public:
  ParseError(const std::string &reason, size_t line)
      : reason(fmt::format("{}: {}", line, reason)) {}
  ParseError(const ParseError &) = default;
  ParseError(ParseError &&) = default;
  const char *what() const noexcept override { return reason.c_str(); }
};

namespace {

struct GraphMLAttribute {
  std::string type;
  std::string name;
};

bool parseBool(const std::string &text) {
  std::string boolString = std::string(text);
  boost::algorithm::to_lower(boolString);
  if (boolString == "false" || boolString == "0" || boolString == "" ||
      boolString == "no")
    return false;
  return true;
}

} // namespace

PowerGrid readGraphML(const std::string &filename, bool all_zero_injection) {

  using namespace std::string_literals;

  PowerGrid outgraph;
  tinyxml2::XMLDocument doc;

  // Load file
  auto err = doc.LoadFile(filename.c_str());
  if (err != tinyxml2::XMLError::XML_SUCCESS) {
    throw ParseError(doc.ErrorIDToName(err), 0);
  }

  // Find graphml
  auto graphml = doc.FirstChildElement("graphml");
  if (!graphml)
    throw ParseError("no graphml content", doc.GetLineNum());

  // Find graph
  auto graph = graphml->FirstChildElement("graph");
  if (!graph)
    throw ParseError("no graph", graphml->GetLineNum());

  // Find attributes
  std::map<std::string, GraphMLAttribute> nodeAttributes;
  for (auto attributes = graphml->FirstChildElement("key");
       attributes != nullptr;
       attributes = attributes->NextSiblingElement("key")) {
    const char *id, *type, *element, *name;
    auto queryAttribute = [&doc, &attributes](const char *name,
                                              const char **attr) -> bool {
      auto err = attributes->QueryAttribute(name, attr);
      if (err) {
        fmt::print(stderr, "[WARN] cannot read attribute {}: {}", name,
                   doc.ErrorIDToName(err));
        return false;
      } else {
        return true;
      }
    };
    if (queryAttribute("id", &id) && queryAttribute("for", &element) &&
        queryAttribute("attr.type", &type) &&
        queryAttribute("attr.name", &name)) {
      nodeAttributes[id] = {type, name};
    }
  }

  // Read vertices
  std::map<std::string, PowerGrid::vertex_descriptor> vertices;
  for (auto node = graph->FirstChildElement("node"); node != nullptr;
       node = node->NextSiblingElement("node")) {
    const char *key;
    bool zero_injection = all_zero_injection;
    if (node->QueryAttribute("id", &key))
      throw ParseError("invalid node", node->GetLineNum());
    std::string name(key), name_alt(key);

    // Read vertex attributes
    for (auto data = node->FirstChildElement("data"); data != nullptr;
         data = data->NextSiblingElement("data")) {
      auto text = data->GetText();
      if (!text)
        throw ParseError("could not read text", data->GetLineNum());
      if (!data->QueryAttribute("key", &key)) {
        if ("zero_injection"s == nodeAttributes[key].name) {
          zero_injection = all_zero_injection || parseBool(text);
        } else if ("name"s == nodeAttributes[key].name) {
          name_alt = std::string(text);
        }
      }
    }

    // Add vertex
    auto vertex =
        boost::add_vertex(Bus{.name = name_alt,
                              .id = static_cast<long>(vertices.size()),
                              .zero_injection = zero_injection},
                          outgraph);
    vertices[name] = vertex;
  }

  // Read edges
  for (auto edge = graph->FirstChildElement("edge"); edge != nullptr;
       edge = edge->NextSiblingElement("edge")) {
    const char *source, *target;
    if (edge->QueryAttribute("source", &source))
      throw ParseError("cannot parse edge source", edge->GetLineNum());
    if (edge->QueryAttribute("target", &target))
      throw ParseError("cannot parse edge target", edge->GetLineNum());
    if (!vertices.contains(source) || !vertices.contains(target))
      throw ParseError("invalid edge", edge->GetLineNum());
    boost::add_edge(vertices[source], vertices[target], outgraph);
  }

  return outgraph;
}

} // end of namespace pds