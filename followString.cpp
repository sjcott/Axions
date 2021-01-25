

// Saving this code in case it becomes useful in the future. This method attempts to connect the strings up and follow it's path through

// Now need to determine which of these coordinates are next to each other. Then can approximate the length by summing these distances
            // For each coordinate try to find coordinates in faces that exist on the same cube. Choose one of the two cubes to check as well.

            double stringLength = 0;
            int fInd = 0; // The index of the point currently following
            for(i=0;i<xString.size();i++){

                bool NextIntersectFound = false;
                double minDistSqr = pow(dx,2) + pow(dy,2) + pow(dz,2); // Start with largest possible distance
                int indClosest;

                for(j=0;j<xString.size();j++){
                    if(j!=fInd){

                        // Check the cube on the +ve side of the face
                        // Might be useful to simplify this with a recursive function?

                        double distSqr = 0;

                        // Check if inside x range

                        if(xString[j] >= dx*floor(xString[fInd]/dx) && xString[j] <= dx*floor(xString[fInd]/dx) + dx){

                            distSqr = pow(xString[j] - xString[fInd],2);

                            // Check if inside y range

                            if(yString[j] >= dy*floor(yString[fInd]/dy) && yString[j] <= dy*floor(yString[fInd]/dy) + dy){

                                distSqr += pow(yString[j] - yString[fInd],2);

                                // Check if inside z range

                                if(zString[j] >= dz*floor(zString[fInd]/dz) && zString[j] <= dz*floor(zString[fInd]/dy) + dy){

                                    distSqr += pow(zString[j] - zString[fInd],2);

                                    NextIntersectFound = true; // If reached this point, this point lies on a face of the next cube
                                    if(distSqr < minDistSqr){ minDistSqr = distSqr;   indClosest = j; }

                                } else if(zString[j] + nz*dz >= dz*floor(zString[fInd]/dz) && zString[fInd] + nz*dz <= dz*floor(zString[fInd]/dz) + dz){

                                    distSqr += pow(zString[j] + nz*dz - zString[fInd],2);

                                    NextIntersectFound = true;
                                    if(distSqr < minDistSqr){ minDistSqr = distSqr;   indClosest = j; }

                                }



                            } else if(yString[j] + ny*dy >= dy*floor(yString[fInd]/dy) && yString[j] + ny*dy <= dy*floor(yString[fInd]/dy) + dy){

                                distSqr += pow(yString[j] + ny*dy - yString[fInd],2);


                                if(zString[j] >= dz*floor(zString[fInd]/dz) && zString[j] <= dz*floor(zString[fInd]/dy) + dy){

                                    distSqr += pow(zString[j] - zString[fInd],2);

                                    NextIntersectFound = true;
                                    if(distSqr < minDistSqr){ minDistSqr = distSqr;   indClosest = j; }

                                } else if(zString[j] + nz*dz >= dz*floor(zString[fInd]/dz) && zString[fInd] + nz*dz <= dz*floor(zString[fInd]/dz) + dz){

                                    distSqr += pow(zString[j] + nz*dz - zString[fInd],2);

                                    NextIntersectFound = true;
                                    if(distSqr < minDistSqr){ minDistSqr = distSqr;   indClosest = j; }

                                }



                            }


                        } else if(xString[j] + nx*dx >= dx*floor(xString[fInd]/dx) && xString[j] + nx*dx <= dx*floor(xString[fInd]/dx) + dx){

                            distSqr = pow(xString[j] + nx*dx - xString[fInd],2);


                            if(yString[j] >= dy*floor(yString[fInd]/dy) && yString[j] <= dy*floor(yString[fInd]/dy) + dy){

                                distSqr += pow(yString[j] - yString[fInd],2);


                                if(zString[j] >= dz*floor(zString[fInd]/dz) && zString[j] <= dz*floor(zString[fInd]/dy) + dy){

                                    distSqr += pow(zString[j] - zString[fInd],2);

                                    NextIntersectFound = true;
                                    if(distSqr < minDistSqr){ minDistSqr = distSqr;   indClosest = j; }

                                } else if(zString[j] + nz*dz >= dz*floor(zString[fInd]/dz) && zString[fInd] + nz*dz <= dz*floor(zString[fInd]/dz) + dz){

                                    distSqr += pow(zString[j] + nz*dz - zString[fInd],2);

                                    NextIntersectFound = true;
                                    if(distSqr < minDistSqr){ minDistSqr = distSqr;   indClosest = j; }

                                }



                            } else if(yString[j] + ny*dy >= dy*floor(yString[fInd]/dy) && yString[j] + ny*dy <= dy*floor(yString[fInd]/dy) + dy){

                                distSqr += pow(yString[j] + ny*dy - yString[fInd],2);


                                if(zString[j] >= dz*floor(zString[fInd]/dz) && zString[j] <= dz*floor(zString[fInd]/dy) + dy){

                                    distSqr += pow(zString[j] - zString[fInd],2);

                                    NextIntersectFound = true;
                                    if(distSqr < minDistSqr){ minDistSqr = distSqr;   indClosest = j; }

                                } else if(zString[j] + nz*dz >= dz*floor(zString[fInd]/dz) && zString[fInd] + nz*dz <= dz*floor(zString[fInd]/dz) + dz){

                                    distSqr += pow(zString[j] + nz*dz - zString[fInd],2);

                                    NextIntersectFound = true;
                                    if(distSqr < minDistSqr){ minDistSqr = distSqr;   indClosest = j; }

                                }



                            }


                        }

                    }
                }

                // Now checked all other intersection points. If NextIntersectFound == false, need to check the -ve side of the face
                if(!NextIntersectFound){
                    for(j=0;j<xString.size();j++){
                        if(j!=fInd){

                            double distSqr = 0;

                            if(xString[j] >= dx*ceil(xString[fInd]/dx) - dx && xString[j] <= dx*ceil(xString[fInd]/dx)){

                                distSqr = pow(xString[j] - xString[fInd],2);


                                if(yString[j] >= dy*ceil(yString[fInd]/dy) - dy && yString[j] <= dy*ceil(yString[fInd]/dy)){

                                    distSqr += pow(yString[j] - yString[fInd],2);


                                    if(zString[j] >= dz*ceil(zString[fInd]/dz) - dz && zString[j] <= dz*ceil(zString[fInd]/dy)){

                                        distSqr += pow(zString[j] - zString[fInd],2);

                                        NextIntersectFound = true; // If reached this point, this point lies on a face of the next cube
                                        if(distSqr < minDistSqr){ minDistSqr = distSqr;   indClosest = j; }

                                } else if(zString[j] + nz*dz >= dz*floor(zString[fInd]/dz) && zString[fInd] + nz*dz <= dz*floor(zString[fInd]/dz) + dz){

                                    distSqr += pow(zString[j] + nz*dz - zString[fInd],2);

                                    NextIntersectFound = true;
                                    if(distSqr < minDistSqr){ minDistSqr = distSqr;   indClosest = j; }

                                }



                            } else if(yString[j] + ny*dy >= dy*floor(yString[fInd]/dy) && yString[j] + ny*dy <= dy*floor(yString[fInd]/dy) + dy){

                                distSqr += pow(yString[j] + ny*dy - yString[fInd],2);


                                if(zString[j] >= dz*floor(zString[fInd]/dz) && zString[j] <= dz*floor(zString[fInd]/dy) + dy){

                                    distSqr += pow(zString[j] - zString[fInd],2);

                                    NextIntersectFound = true;
                                    if(distSqr < minDistSqr){ minDistSqr = distSqr;   indClosest = j; }

                                } else if(zString[j] + nz*dz >= dz*floor(zString[fInd]/dz) && zString[fInd] + nz*dz <= dz*floor(zString[fInd]/dz) + dz){

                                    distSqr += pow(zString[j] + nz*dz - zString[fInd],2);

                                    NextIntersectFound = true;
                                    if(distSqr < minDistSqr){ minDistSqr = distSqr;   indClosest = j; }

                                }



                            }


                        } else if(xString[j] + nx*dx >= dx*floor(xString[fInd]/dx) && xString[j] + nx*dx <= dx*floor(xString[fInd]/dx) + dx){

                            distSqr = pow(xString[j] + nx*dx - xString[fInd],2);


                            if(yString[j] >= dy*floor(yString[fInd]/dy) && yString[j] <= dy*floor(yString[fInd]/dy) + dy){

                                distSqr += pow(yString[j] - yString[fInd],2);


                                if(zString[j] >= dz*floor(zString[fInd]/dz) && zString[j] <= dz*floor(zString[fInd]/dy) + dy){

                                    distSqr += pow(zString[j] - zString[fInd],2);

                                    NextIntersectFound = true;
                                    if(distSqr < minDistSqr){ minDistSqr = distSqr;   indClosest = j; }

                                } else if(zString[j] + nz*dz >= dz*floor(zString[fInd]/dz) && zString[fInd] + nz*dz <= dz*floor(zString[fInd]/dz) + dz){

                                    distSqr += pow(zString[j] + nz*dz - zString[fInd],2);

                                    NextIntersectFound = true;
                                    if(distSqr < minDistSqr){ minDistSqr = distSqr;   indClosest = j; }

                                }



                            } else if(yString[j] + ny*dy >= dy*floor(yString[fInd]/dy) && yString[j] + ny*dy <= dy*floor(yString[fInd]/dy) + dy){

                                distSqr += pow(yString[j] + ny*dy - yString[fInd],2);


                                if(zString[j] >= dz*floor(zString[fInd]/dz) && zString[j] <= dz*floor(zString[fInd]/dy) + dy){

                                    distSqr += pow(zString[j] - zString[fInd],2);

                                    NextIntersectFound = true;
                                    if(distSqr < minDistSqr){ minDistSqr = distSqr;   indClosest = j; }

                                } else if(zString[j] + nz*dz >= dz*floor(zString[fInd]/dz) && zString[fInd] + nz*dz <= dz*floor(zString[fInd]/dz) + dz){

                                    distSqr += pow(zString[j] + nz*dz - zString[fInd],2);

                                    NextIntersectFound = true;
                                    if(distSqr < minDistSqr){ minDistSqr = distSqr;   indClosest = j; }

                                }



                            }


                        }

                    }
                }


            }